! Copyright (c) 2019-2026 Armin Corbin
! University of Bonn, Institute for Geodesy and Geoinformation,
! Astronomical, Physical and Mathematical Geodesy (APMG) group
!
! This file is part of pdaf-binding-tiegcm, the coupling layer between
! the Parallel Data Assimilation Framework (PDAF) and TIE-GCM.
!
! pdaf-binding-tiegcm is free software: you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public License
! as published by the Free Software Foundation, either version 3 of
! the License, or (at your option) any later version.
!
! pdaf-binding-tiegcm is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with pdaf-binding-tiegcm. If not, see
! <http://www.gnu.org/licenses/>.
!
! This software is part of the NCAR TIE-GCM. Use is governed by the Open
! Source Academic Research License Agreement contained in the file
! tiegcmlicense.txt.
!
! Reads along-track satellite observation/trajectory files (TOLEOS, GROOPS, DENSWIND) and spline-interpolates them for use as PDAF observations.

module trajectory_data_module

  use interpolation_module, only: spline_interpolator

  implicit none

  type :: trajectory_data
    ! one could reuse the jacobian. But the initailization is fast anyway so this has not a high priority
    real, allocatable, dimension(:) :: modeltime ! seconds since model reference epoch
    logical, dimension(:), allocatable :: data_gap
    type(spline_interpolator) :: values
    type(spline_interpolator) :: x_trf
    type(spline_interpolator) :: y_trf
    type(spline_interpolator) :: z_trf
    real, dimension(3) :: current_position
    real :: current_value
    real, private :: time_previous_call = -huge(0.0)
    contains
    procedure, pass(this) :: deallocate => trajectory_data_deallocate
    procedure, pass(this) :: read_denswind => trajectory_data_read_denswind
    procedure, pass(this) :: read_txt => trajectory_data_read_txt
    procedure, pass(this) :: get_at_epoch => trajectory_data_get_at_epoch
    procedure, pass(this) :: valid_at_epoch => trajectory_data_valid_at_epoch
  end type

  contains

  !> Deallocates a trajectory's time/data-gap arrays and its spline interpolators.
  subroutine trajectory_data_deallocate(this)

    implicit none

    ! arguments
    class(trajectory_data) :: this

    if(allocated(this%modeltime)) deallocate(this%modeltime)
    if(allocated(this%data_gap)) deallocate(this%data_gap)

    call this%values%deallocate
    call this%x_trf%deallocate
    call this%y_trf%deallocate
    call this%z_trf%deallocate

  end subroutine

  !> Reads one or more along-track text files (TOLEOS or GROOPS format), converts
  !! positions to cartesian coordinates, and builds spline interpolators plus a
  !! data-gap mask.
  subroutine trajectory_data_read_txt(this, file, file_format)

  ! intern
  use coordinates_module, only: wgs84
  use time_module, only: convert_GPS_to_UTC, TIME_SYSTEM_GPS, convert_to_seconds_since_ref
  use configuration, only: path_len, cfg_log

  ! extern
  use ESMF

  implicit none

  ! arguments
  class(trajectory_data) :: this
  character(len=path_len), intent(in) :: file
  character(len=*), intent(in) :: file_format

  ! local
  character(len=path_len), dimension(:), allocatable :: matches
  integer :: nfiles

  integer :: i
  integer :: n_total_lines
  integer :: n_valid_lines
  integer :: n_obs

  character(len=27) :: time_string


  type(ESMF_Time), dimension(:), allocatable :: time
  integer, dimension(:), allocatable :: time_system
  real, dimension(:,:), allocatable :: data

  real, allocatable, dimension(:,:) :: time_diff

  write(*,*) '------------------------------------------------------------------------------'
  write(*,*) 'read along track data from ', file, " expecting the format: ", file_format

  ! ATTENTION we assume that alphabetic and chronological order are the same
  call get_matching_files(file, matches)

  nfiles = size(matches,dim=1)

  n_total_lines = 0
  do i = 1, nfiles
    n_total_lines = n_total_lines + count_lines(matches(i))
  end do

  allocate(time(n_total_lines))
  allocate(time_system(n_total_lines))
  allocate(data(n_total_lines,4))

  n_obs = 0
  do i = 1, nfiles
    select case(trim(file_format))
    case('toleos_reduced','toleos','toleos_short')
      call parse_toleos(matches(i),file_format,time,time_system,data,n_valid_lines,offset=n_obs)
    case('groops')
      call parse_groops(matches(i),time,time_system,data,n_valid_lines,offset=n_obs)
    case default
      call shutdown(trim(file_format)//" is an invalid toleos format")
    end select

    n_obs = n_obs + n_valid_lines
    if(cfg_log%verbose_level>0) then
      write(*,*) "file: ", i, ". Number of observations in this file: ",&
                 n_valid_lines, " total number so far: ",  n_obs
      write(*,*) ".............................................................................."
    end if
  end do

  deallocate(matches)

  call convert_GPS_to_UTC(time(1:n_obs), time_system(1:n_obs)==TIME_SYSTEM_GPS)

  allocate(this%modeltime(n_obs))
  call convert_to_seconds_since_ref(time(1:n_obs),this%modeltime)

  ! convert to cartesian coordinates to avoid jumps in longitude interpolation
  do i=1,n_obs
      data(i,1:3) = wgs84%ellipsoidalGeodetic2cartesian(data(i,1:3))
  end do

  ! initalize spline interpolation
  call this%x_trf%init(this%modeltime,data(1:n_obs,1), degree=3)
  call this%y_trf%init(this%modeltime,data(1:n_obs,2), degree=3)
  call this%z_trf%init(this%modeltime,data(1:n_obs,3), degree=3)
  call this%values%init(this%modeltime,data(1:n_obs,4), degree=3)

  ! check for data gaps in time series
  allocate(this%data_gap(n_obs))
  this%data_gap = .false.

  ! mark observations adjacent to data gaps with flag
  allocate(time_diff(2,n_obs))
  time_diff = 0.0

  time_diff(1,1:n_obs-1)=this%modeltime(2:n_obs)-this%modeltime(1:n_obs-1)
  time_diff(2,2:n_obs)=time_diff(1,1:n_obs-1)

  ! consider everything longer than 30s as data gap
  where( maxval(time_diff,dim=1)>30.0) this%data_gap = .true.

  ! extend gap to neigbour on left and right twice
  ! so bspline of degree is always in valid range
  call binary_dilation(this%data_gap,2)

  deallocate(time_diff)

  if(cfg_log%verbose_level>0) then
    write(*,*) "print the first 10 parsed and processed observations:"
    write(*,'(a21, a4, 4a15)') "time","gap","x (m)","y (m)","z (m)","den (kg m-3)"
    do i = 1, min(10, n_obs)
      call ESMF_TimeGet(time(i), timeString=time_string)
      write(*,'(a21, l3, 3f15.1, e15.8)')  time_string, this%data_gap(i), data(i,:)
    end do
  end if

  deallocate(time)
  deallocate(time_system)
  deallocate(data)

  end subroutine

  !> Counts the number of lines in a text file.
  function count_lines(file) result(n_lines)

    implicit none

    ! arguments
    character(len=*), intent(in) :: file

    ! result
    integer :: n_lines

    ! local
    integer :: io
    integer :: rc

    open(newunit=io, file=file, status="old", action="read")
    n_lines = 0
    do
      read (io, *,iostat=rc)
      if (rc /= 0) exit
      n_lines = n_lines+1
    end do
    close(io)
  end function

  !> Parses one TOLEOS-format along-track file into time/position/value arrays,
  !! optionally skipping flagged anomalous records.
  subroutine parse_toleos(file,toleos_format,time,time_system,data,n_valid_lines,&
                          offset,&
                          ignore_anomalous_data)

    use time_module, only: iso_time_str_to_esfm_time, TIME_SYSTEM_UTC, TIME_SYSTEM_GPS
    use configuration, only: cfg_log

    ! extern
    use ESMF

    implicit none

    ! arguments
    character(len=*), intent(in) :: file
    character(len=*), intent(in) :: toleos_format
    type(ESMF_Time), dimension(:), intent(inout) :: time
    integer, dimension(:), intent(out) :: time_system
    real, dimension(:,:), intent(inout) :: data ! lon, lat, alt, density
    integer, intent(out) :: n_valid_lines
    integer, intent(in), optional :: offset
    logical, intent(in), optional :: ignore_anomalous_data

    ! local
    integer :: offset_
    logical :: ignore_anomalous_data_

    logical :: exists
    integer :: io
    integer :: rc
    integer :: i
    integer :: anomalous_data_count
    logical :: has_flag

    character(len=27) :: time_string

    character(len=128) :: line_format

    real :: flag

    if(present(offset))then
      offset_ = offset
    else
      offset_ = 0
    end if

    if(present(ignore_anomalous_data))then
      ignore_anomalous_data_ = ignore_anomalous_data
    else
      ignore_anomalous_data_ = .true.
    end if

    anomalous_data_count = 0

    inquire(file=file, exist=exists)

    select case(trim(toleos_format))
    case('toleos_reduced')
      line_format ='(a27,1x,f10.3,1X, f8.3,1X, f7.3,1X,  6X,1X,   7X,1X,e15.8,1X,  15X,1X,f4.1,1X,f4.1)'
      has_flag = .true.
    case('toleos')
      line_format ='(a27,1x,f10.3,1X,f13.8,1X,f13.8,1X,  6X,1X,  13X,1X,e15.8,1X,  15X,1X,f4.1,1X,4X)'
      has_flag = .true.
    case('toleos_short')
      line_format ='(a27,1x,f10.3,1X,f13.8,1X,f13.8,1X,  6X,1X,  13X,1X,e15.8)'
      has_flag = .false.
    case default
      call shutdown(trim(toleos_format)//" is an invalid toleos format")
      has_flag = .false.
    end select

    if (exists) then
      open(newunit=io, file=file, status="old", action="read",form='FORMATTED')
      i = offset_+1
      do
        if(has_flag) then
          read (io, line_format, &
                iostat=rc) time_string, data(i,3), data(i,1), data(i,2), data(i,4), flag
        else
          read (io, line_format, &
              iostat=rc) time_string, data(i,3), data(i,1), data(i,2), data(i,4)
        end if
        if (rc < 0) exit
        if ((time_string(1:1).eq.'#').or.(len_trim(time_string)==0)) cycle

        if((has_flag.eqv..true.) .and. (flag /= 0)) then
          anomalous_data_count = anomalous_data_count + 1
          if(ignore_anomalous_data_) cycle
        end if

        call iso_time_str_to_esfm_time(time_string(1:21),time(i))

        select case(time_string(25:27))
          case("UTC")
            time_system(i)=TIME_SYSTEM_UTC
          case("GPS")
            time_system(i)=TIME_SYSTEM_GPS
          case default
            write(*,*) time_string(24:27), " is not a supported time system"
        end select

!          call ESMF_TimeGet(time(i), timeString=time_string)
!          write(*,*) i, time_string, data(i,:)

        i = i+1

      end do
      n_valid_lines = i-offset_-1

      write(*,*) 'parsed ', n_valid_lines, ' lines', &
                  merge(' ignored  ',' included ',ignore_anomalous_data_), &
                  anomalous_data_count, ' anomalous data points'

      close (io)

      if(cfg_log%verbose_level>0) then
        write(*,*) "print the first 10 parsed observations:"
        write(*,'(a21,3a12,a15)') "time","lon (deg)","lat (deg)","alt (m)","den (kg m-3)"
        do i = offset_+1, min(offset_+11, offset_+n_valid_lines)
          call ESMF_TimeGet(time(i), timeString=time_string)
          write(*,'(a21, 2f12.6, f12.1, e15.8)') time_string,  data(i,:)
        end do
      end if
    else
      write(*,*) trim(file), 'does not exist'
      n_valid_lines = 0
    end if
  end subroutine

  !> Parses one GROOPS-format along-track file (modified-Julian-date based) into
  !! time/position/value arrays.
  subroutine parse_groops(file,time,time_system,data,n_valid_lines,&
                          offset)

    use time_module, only: mjd_to_esmf_time, TIME_SYSTEM_GPS
    use configuration, only: cfg_log

    ! extern
    use ESMF

    implicit none

    ! arguments
    character(len=*), intent(in) :: file
    type(ESMF_Time), dimension(:), intent(inout) :: time
    integer, dimension(:), intent(out) :: time_system
    real, dimension(:,:), intent(inout) :: data ! lon, lat, alt, density
    integer, intent(out) :: n_valid_lines
    integer, intent(in), optional :: offset


    ! local
    integer :: offset_

    logical :: exists
    integer :: io
    integer :: rc
    integer :: i

    real :: mjd

    character(len=128) :: line_format
    character(len=27) :: time_string


    if(present(offset))then
      offset_ = offset
    else
      offset_ = 0
    end if

    inquire(file=file, exist=exists)
    line_format =  '(F25.18,1X,F25.18,1X,F25.18,1X,F25.18,1X,25X,1X,F25.18)'

    time_system=TIME_SYSTEM_GPS

    if (exists) then
      open(newunit=io, file=file, status="old", action="read",form='FORMATTED')

      do i = 1, 6
          read(io, *, iostat=rc)
          if (rc /= 0) error stop 'File ended while skipping header lines'
      end do

      i = offset_+1
      do
        read (io, line_format, &
              iostat=rc) mjd, data(i,1), data(i,2), data(i,3), data(i,4)
        if (rc < 0) exit

        time(i) = mjd_to_esmf_time(mjd)

        i = i+1

      end do
      n_valid_lines = i-offset_-1

      write(*,*) 'parsed ', n_valid_lines, ' lines'

      close (io)

      if(cfg_log%verbose_level>0) then
        write(*,*) "print the first 10 parsed observations:"
        write(*,'(a21,3a12,a15)') "time","lon (deg)","lat (deg)","alt (m)","den (kg m-3)"
        do i = offset_+1, min(offset_+11, offset_+n_valid_lines)
          call ESMF_TimeGet(time(i), timeString=time_string)
          write(*,'(a21, 2f12.6, f12.1, e15.8)') time_string,  data(i,:)
        end do
      end if
    else
      write(*,*) trim(file), 'does not exist'
      n_valid_lines = 0
    end if
  end subroutine

  !> Reads one or more DENSWIND NetCDF files, converts positions to cartesian
  !! coordinates, and builds spline interpolators plus a data-gap mask.
  subroutine trajectory_data_read_denswind(this, file)

    ! extern
    use netcdf
    use ESMF

    ! tie-gcm
    use nchist_module, only: handle_ncerr

    ! intern
    use coordinates_module, only: wgs84
    use netcdf_functionality, only: read_time_variable
    use time_module, only: convert_GPS_to_UTC, convert_to_seconds_since_ref

    implicit none

    ! arguments
    class(trajectory_data) :: this
    character(len=*), intent(in) :: file

    ! local
    integer :: istat
    integer, dimension(:), allocatable :: ncid

    integer :: var_id

    type(ESMF_TIME), allocatable, dimension(:) :: time, time_all

!     character(len=64) :: timestr

    real, allocatable, dimension(:) :: values
    real, allocatable, dimension(:,:) :: position

    integer, dimension(:), allocatable :: counts
    integer, dimension(:), allocatable :: lbounds
    integer, dimension(:), allocatable :: ubounds
    integer :: n_all

    integer :: i

    integer :: nfiles
    integer :: fid
    integer :: dim_id

    character(len=256), dimension(:), allocatable :: matches

    write(*,*) 'read along track data from ', file

    ! ATTENTION we assume that alphabetic and chronological order are the same
    call get_matching_files(file, matches)

    nfiles = size(matches,dim=1)

    allocate(ncid(nfiles))
    allocate(counts(nfiles))
    allocate(lbounds(nfiles))
    allocate(ubounds(nfiles))

    do fid = 1, nfiles
      istat = nf90_open(path=trim(matches(fid)), &
                        mode=NF90_NOWRITE, &
                        ncid=ncid(fid))
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'cannot open '//matches(fid))

      istat = nf90_inq_dimid( ncid(fid), 'time', dim_id )
      istat = nf90_inquire_dimension( ncid(fid), dim_id, len=counts(fid) )
    end do

    ! determine position of each file in array containing data from all files
    n_all = sum(counts)
    lbounds(1) = 1
    do fid = 2, nfiles
      lbounds(fid) = lbounds(fid-1) + counts(fid-1)
    end do
    ubounds = lbounds + counts -1

    ! read time
    allocate(time_all(n_all))

    do fid = 1, nfiles
      call read_time_variable(ncid(fid), time)
      time_all(lbounds(fid):ubounds(fid)) = time
      deallocate(time)
    end do

    call convert_GPS_to_UTC(time_all)
    allocate(this%modeltime(n_all))
    call convert_to_seconds_since_ref(time_all,this%modeltime)

    ! read data gaps
    allocate(this%data_gap(n_all))
    do fid = 1, nfiles
      call denswind_get_data_gap(ncid(fid), this%data_gap(lbounds(fid):ubounds(fid)))
    end do

    ! read trajectory
    allocate(values(n_all))
    allocate(position(n_all,3))

    do fid = 1, nfiles
      istat = nf90_inq_varid(ncid=ncid(fid), name="density", varid=var_id)
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'trajectory_data_read_denswind')
      istat = nf90_get_var(ncid=ncid(fid), varid=var_id, values=values(lbounds(fid):ubounds(fid)))
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'trajectory_data_read_denswind')

      istat = nf90_inq_varid(ncid=ncid(fid), name="lon", varid=var_id)
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'trajectory_data_read_denswind')
      istat = nf90_get_var(ncid=ncid(fid), varid=var_id, values=position(lbounds(fid):ubounds(fid),1))
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'trajectory_data_read_denswind')

      istat = nf90_inq_varid(ncid=ncid(fid), name="lat", varid=var_id)
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'trajectory_data_read_denswind')
      istat = nf90_get_var(ncid=ncid(fid), varid=var_id, values=position(lbounds(fid):ubounds(fid),2))
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'trajectory_data_read_denswind')

      istat = nf90_inq_varid(ncid=ncid(fid), name="alt", varid=var_id)
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'trajectory_data_read_denswind')
      istat = nf90_get_var(ncid=ncid(fid), varid=var_id, values=position(lbounds(fid):ubounds(fid),3))
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'trajectory_data_read_denswind')
    end do

    this%data_gap = (this%data_gap .or. &
                     isnan(values) .or. &
                     isnan(position(:,1)) .or. &
                     isnan(position(:,2)) .or. &
                     isnan(position(:,3)))

!     do i=1,n_all
!       if(this%data_gap(i).eqv..false.)then
!         call ESMF_TimeGet(time_all(i), timeString=timestr)
!         write(*,*) timestr, position(i,:), values(i)
!       end if
!     end do

    deallocate(time_all)

    ! convert to cartesian coordinates to avoid jumps in longitude interpolation
    do i=1,n_all
      if(this%data_gap(i).eqv..false.)then
        position(i,:) = wgs84%ellipsoidalGeodetic2cartesian(position(i,:))
      end if
    end do

    ! initalize spline interpolation
    call this%x_trf%init( pack(this%modeltime,this%data_gap.eqv..false.),&
                          pack(position(:,1),this%data_gap.eqv..false.),&
                          degree=3)
    call this%y_trf%init( pack(this%modeltime,this%data_gap.eqv..false.),&
                          pack(position(:,2),this%data_gap.eqv..false.),&
                          degree=3)
    call this%z_trf%init( pack(this%modeltime,this%data_gap.eqv..false.),&
                          pack(position(:,3),this%data_gap.eqv..false.),&
                          degree=3)
    call this%values%init( pack(this%modeltime,this%data_gap.eqv..false.),&
                          pack(values,this%data_gap.eqv..false.),&
                          degree=3)

    deallocate(values)
    deallocate(position)

    ! extend gap to neigbour on left and right twice
    ! so bspline of degree is always in valid range
    call binary_dilation(this%data_gap,2)

    deallocate(ncid,counts,lbounds,ubounds)

  end subroutine

  !> Reads the accelerometer/orbit/star-camera gap flags from a DENSWIND file and
  !! combines them into a single data-gap mask.
  subroutine denswind_get_data_gap(ncid,data_gap)

    ! extern
    use netcdf

    ! tie-gcm
    use nchist_module, only: handle_ncerr

    implicit none

    ! arguments
    integer, intent(in) :: ncid
    logical, dimension(:), intent(inout) :: data_gap

    ! local
    integer, allocatable, dimension(:,:) :: gaps

    integer :: gap_grp_id

    integer :: istat, var_id
    integer :: n

    n = size(data_gap)

    allocate(gaps(n,3))

    istat = nf90_inq_ncid(ncid, 'gap', gap_grp_id)
    istat = nf90_inq_varid(ncid=gap_grp_id, name="accelerometer", varid=var_id)
    istat = nf90_get_var(ncid=gap_grp_id, varid=var_id, values=gaps(:,1))
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'cannot read gaps in accelerometer')

    istat = nf90_inq_varid(ncid=gap_grp_id, name="orbit", varid=var_id)
    istat = nf90_get_var(ncid=gap_grp_id, varid=var_id, values=gaps(:,2))
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'cannot read gaps in orbit')

    istat = nf90_inq_varid(ncid=gap_grp_id, name="starcamera", varid=var_id)
    istat = nf90_get_var(ncid=gap_grp_id, varid=var_id, values=gaps(:,3))
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'cannot read gaps in starcamera')

    data_gap = sum(gaps,dim=2) /= 0

    deallocate(gaps)

  end subroutine

  !> Grows a logical mask by the given number of iterations, marking neighbours of
  !! true values as true (used to widen data-gap regions).
  subroutine binary_dilation(mask,iterations)
    implicit none

    ! arguments
    logical, dimension(:), intent(inout) :: mask
    integer, intent(in) :: iterations

    ! local
    integer :: n
    integer :: i

    n = size(mask)

    do i=1,iterations
      mask(1:n-1) = mask(1:n-1) .or. mask(2:n)
      mask(2:n) = mask(2:n) .or. mask(1:n-1)
    end do
  end subroutine

  !> Returns whether a trajectory has valid (non-data-gap, in-range) observations
  !! at the given model time.
  function trajectory_data_valid_at_epoch(this,modeltime) result(res)

    use search_module, only: nearest_neighbour, NN_RIGHT

    implicit none

    ! arguments
    class(trajectory_data) :: this
    real, intent(in) :: modeltime

    ! result
    logical :: res

    ! local
    integer :: idx

    res = .true.

    if(modeltime < this%modeltime(1)) then
      write(*,*) 'requested time is befor interval of observation'
      res = .false.
      return
    end if

    if(modeltime > this%modeltime(size(this%modeltime,dim=1))) then
      write(*,*) 'requested time is after interval of observation'
      res = .false.
      return
    end if

    idx = nearest_neighbour(this%modeltime,modeltime,NN_RIGHT)
    if(this%data_gap(idx)) res = .false.
    if(idx>1) then
      if(this%data_gap(idx-1)) res = .false.
    end if
  end function

  !> Spline-interpolates and caches the trajectory's value and geodetic position at
  !! the given model time.
  subroutine trajectory_data_get_at_epoch(this,modeltime,val,position)

    use angle_module, only: map_angle_to_interval_minus_180_and_180
    use coordinates_module, only: wgs84

    implicit none

    ! arguments
    class(trajectory_data) :: this
    real, intent(in) :: modeltime
    real, intent(out) :: val
    real, dimension(3), intent(out) :: position

    if(this%time_previous_call /= modeltime) then
      ! TODO high multiplicity knots at gaps?

      this%current_value = this%values%interpolate(modeltime)
      this%current_position(1) = this%x_trf%interpolate(modeltime)
      this%current_position(2) = this%y_trf%interpolate(modeltime)
      this%current_position(3) = this%z_trf%interpolate(modeltime)

      ! TODO which ellipsoid is used in TIE-GCM
      this%current_position = wgs84%cartesian2ellipsoidalGeodetic(this%current_position)
      this%current_position = wgs84%geodetic2geocentricLatitude(this%current_position)

      this%current_position(1) = map_angle_to_interval_minus_180_and_180(this%current_position(1))

      this%time_previous_call = modeltime
    end if

    val = this%current_value
    position = this%current_position

  end subroutine

  !> Lists files matching a shell wildcard pattern (via `ls`) and returns the
  !! matching filenames, broadcasting the result from MPI rank 0.
  subroutine get_matching_files(file_pattern, matches, comm)

    use mpi_f08

    implicit none

    ! arguments
    character(len=*), intent(in) :: file_pattern ! path to file(s) using sh wildcards, e.g. *
    character(len=256), dimension(:), allocatable, intent(out) :: matches
    type(MPI_Comm), optional, intent(in) :: comm

    ! local
    type(MPI_Comm) :: comm_
    integer :: rank
    integer :: status, io_error, iostat
    integer :: n, i

    character(len=256) :: line
    integer, parameter :: io_unit=20

    if(present(comm))then
      comm_ = comm
    else
      comm_ = MPI_Comm_World
    end if

    call MPI_COMM_RANK(comm_,rank)

    if(rank==0)then
      call system('ls -1 '//file_pattern//' > matching_files.txt', status)
      if(status/=0)then
        write(*,*) 'error listing ', file_pattern
      end if
    end if

    call MPI_Barrier(comm_)
    open(unit=io_unit,&
         file='matching_files.txt',&
         status='old',&
         action='read', &
         iostat = io_error)

    n = 0
    if ( io_error == 0) then

      ! determine line number
      do
        read(unit=io_unit,FMT='(a)',iostat=iostat) line
        if (iostat/=0) EXIT
        n = n+1
      end do

      allocate(matches(n))

      rewind(io_unit)

      do i = 1,n
        read(io_unit,FMT='(a)') matches(i)
      end do
    end if

    close(unit=io_unit)

    write(*,'(a,i4,2a)') 'found ', n, ' files matching the pattern ', trim(file_pattern)
    do i = 1,n
      write(*,'(2x,a)') trim( matches(i))
    end do
  end subroutine

end module
