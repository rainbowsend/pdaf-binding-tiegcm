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
! Debug helpers: synthetic test grids/values plus an instant per-step NetCDF writer for diagnosing crashes.

!> Prints a string/message prefixed with the model instance and rank of the corresponding PE and the current model step
subroutine print_loc(str)
  implicit none

  character(len=*), intent(in) :: str

  character(len=128) :: loc

  call error_source_string(loc)
  write(*,'(a,1x,a)') trim(loc), trim(str)

end subroutine print_loc

!> Formats a string identifying the current ensemble member, world rank, model-internal rank, and time step.
subroutine error_source_string( str )

  ! tie-gcm
  use mpi_module, only: mytid
  use init_module,only: istep

  ! intern
  use mod_parallel_pdaf, only: task_id, rank_world

  implicit none

  character(len=128), intent(out) :: str

  write(str, '(a,i4,1x,a,i4,1x,a,i4,1x,a,i5)') &
                  'model instance: ', task_id, &
                    'world rank:' , rank_world, &
                  'model intern rank: ', mytid, &
                              'step: ', istep
end subroutine error_source_string

!> Prints string prefixed with the current iteration count, step count, and model time.
subroutine print_time_info(str)

  ! tie-gcm
  use init_module, only: istep, iter
  use hist_module,only: modeltime

  implicit none

  character(len=*), intent(in) :: str

  character(len=16), dimension(3) :: formatted

  write( formatted(1), '(i9)'  ) istep
  write( formatted(2), '(i9)'  ) iter

  write(*, '(a,1x,2a,3x,2a,3x,a,i4,3i3)') &
        trim(str), &
        'iteration: ',  trim(adjustl( formatted(1) )) , &
        'steps since 01 JAN 00:00: ' ,  trim(adjustl( formatted(2) )), &
        'model time is', modeltime(1), modeltime(2:4)

end subroutine

module test_values

  implicit none
  contains

  !> Generates a synthetic test value at (lon,lat,h), used to verify interpolation routines.
  function synthetic_val( lon, lat, h ) result( val )

    implicit none

    ! arg
    real, intent(in) :: lon, lat, h

    ! res
    real(kind=8) :: val

    ! local
    real :: x,y,z
    real, parameter :: a = 1
    real, parameter :: b = 2
    real, parameter :: c = 3
    real, parameter :: d2r = 3.14159265358979/180.0

    x = h * a * cos(lon*d2r) * sin(lat*d2r)
    y = h * b * sin(lon*d2r) * sin(lat*d2r)
    z = h * c * cos(lat*d2r)

    val = (x + y + z)*sin(lon*d2r/4)

  end function synthetic_val

  !> Generates a synthetic test value that is linear in lon, lat, and h, used to verify interpolation routines.
  function synthetic_val_lin_fun( lon, lat, h ) result( val )

    implicit none

    ! arg
    real, intent(in) :: lon, lat, h

    ! res
    real(kind=8) :: val

    ! local

    val =3*lon-lat+h/100

  end function synthetic_val_lin_fun

  !> Generates a synthetic density-like test value with day/night and pole/equator variation, used to verify interpolation routines.
  function synthetic_val_density( lon, lat, h ) result( val )

    implicit none

    ! arg
    real, intent(in) :: lon, lat, h

    ! res
    real(kind=8) :: val

    real, parameter :: rho0 = 1E-12 ! g/cm^3
    real, parameter :: h0 = 100E+3 ! m
    real, parameter :: scale_height = 30E+3 ! m
    real, parameter :: day_night_factor = 0.1
    real, parameter :: pol_factor = -0.1;

    real, parameter :: d2r = 3.14159265358979/180.0

    ! local

    val = rho0 * exp(-(h-h0)/scale_height) &
               * (1+day_night_factor*sin(lon*d2r/2)**2) &
               * (1+pol_factor*sin(lat*d2r)**2)


  end function synthetic_val_density

  !> Fills a lon/lat/alt grid with synthetic_val test values.
  function synthetic_grid( lons, lats, alts ) result( syn_grid )

    implicit none

    ! arg
    real, dimension(:), intent(in) :: lons, lats, alts

    ! res
    real(kind=8), dimension( size(lons), size(lats), size(alts) ) :: syn_grid

    ! local
    integer :: i, j, k

    do i=1, size(lons)
        do j=1, size(lats)
            do k=1, size(alts)
                syn_grid(i,j,k) =  synthetic_val( lons(i), lats(j), alts(k)  )
            end do
        end do
    end do

  end function synthetic_grid

end module test_values

module nc_debug

! Tiegcm collects data and writes it at the end of the main loop. Thus, in case of an error data is missing
! This module instantly writes the data so it can be used to debug errors with crashes
! writes specified fields of all ensembler members to the same file. Must be called by all ranks.

USE nchist_module, ONLY: handle_ncerr

use mpi_f08

use mod_parallel_pdaf, only: rank_world

use array_print_module, only: printMat

implicit none

#include <netcdf.inc>

integer :: var_id_time
integer :: var_id_mtime

integer :: istat, ierr
integer :: ncid

integer, parameter :: maxDebugFields = 30
character(len=32), dimension(maxDebugFields) :: var_list
integer, dimension(maxDebugFields) :: var_id
integer :: var_list_len = 0

integer :: ulimlen = 0

character(len=256) :: ncfile 

logical, save :: panic_in_progress = .false.
logical, save :: panic_is_enabled = .false.

contains

!> Writes a warning for every TIE-GCM field that contains NaN values at the current or previous time level.
subroutine nan_check_f4d(description)

  use fields_module, only: itc, itp, f4d

  implicit none

  ! arguments
  character(len=*), intent(in) :: description

  ! local
  integer :: i

  do i = 1, size(f4d)
    if(any(isnan(f4d(i)%data(:,:,:,itc)))) then
      write(*,*) 'nan in ', f4d(i)%short_name, ' itc ', description
    end if
    if(any(isnan(f4d(i)%data(:,:,:,itp)))) then
      write(*,*) 'nan in ', f4d(i)%short_name, ' itp ', description
    end if
  end do

end subroutine

!> Writes every TIE-GCM prognsotic field of the calling rank into a single, self
!! describing NetCDF file, using no MPI communication at all.
!!
!! Intended for use from the shutdown subroutine (util.F), where in general
!! only a subset of the ranks is present
!!
!! The file is named panic_<member>_<rank>_<istep>.nc so that the dumps of
!! different ranks can not overwrite each other.
!!
subroutine panic_write

  ! extern
  ! this module uses the F77 interface (#include <netcdf.inc>, nf_* names);
  ! import the F90 names used here explicitly so the two can not clash
  use netcdf, only: nf90_create, nf90_def_var, nf90_put_att, nf90_enddef, &
                    nf90_put_var, nf90_close, &
                    NF90_NOERR, NF90_NETCDF4, NF90_DOUBLE, NF90_GLOBAL

  ! tie-gcm
  use fields_module, only: f4d, f3d, itc, itp
  use nchist_module, only: handle_ncerr
  use params_module, only: nlevp1
  use init_module, only: istep
  use mpi_module, only: mytid

  ! intern
  use state_module, only: idx_intern, idx_nc
  use mod_parallel_pdaf, only: task_id
  use netcdf_functionality, only: init_spacial_dims_subdomain, &
                                  write_spacial_dims_subdomain, &
                                  add_creation_time

  implicit none

  ! local
  integer :: i, k, istat, ncid, lev_dim
  integer :: nlons, nlats
  integer :: dim_ids(5), coord_var_ids(5)
  integer, allocatable :: var_ids(:)
  integer, allocatable :: var_ids_3d(:)
  integer :: time_slot(2), start_p(4), count_p(4)
  real, dimension(:,:,:), allocatable :: field_reshaped
  character(len=80) :: filename
  character(len=16) :: units

  ! panic write is enalbled after model initalization
  if(.not.panic_is_enabled) return

  ! guards against re-entering while already dumping (shutdown may be reached
  ! again from within this routine)
  if(panic_in_progress) return
  panic_in_progress = .true.

  if(.not.allocated(idx_intern)) then
    write(*,*) 'panic_write: field indices not initialized yet - nothing to dump'
    return
  end if

  nlons = idx_intern(mytid)%nlons
  nlats = idx_intern(mytid)%nlats
  time_slot = (/itp, itc/)

  write(filename,'(a,i0.3,a1,i0.3,a1,i0.5,a3)') 'panic_',task_id,'_',mytid,'_',istep,'.nc'

  istat = nf90_create(filename, NF90_NETCDF4, ncid)
  if (istat /= NF90_NOERR) then
    call handle_ncerr(istat,'panic_write: Error creating '//trim(filename))
    return
  end if

  ! Record where this subdomain sits on the global grid, so that the dumps of
  ! the individual ranks can be interpreted without knowing the decomposition.
  istat = nf90_put_att(ncid, NF90_GLOBAL, 'title', 'TIE-GCM-PDAF panic dump (rank local)')
  istat = nf90_put_att(ncid, NF90_GLOBAL, 'ensemble_member', task_id)
  istat = nf90_put_att(ncid, NF90_GLOBAL, 'rank', mytid)
  istat = nf90_put_att(ncid, NF90_GLOBAL, 'model_step', istep)
  istat = nf90_put_att(ncid, NF90_GLOBAL, 'lon0', idx_nc(mytid)%lon0)
  istat = nf90_put_att(ncid, NF90_GLOBAL, 'lon1', idx_nc(mytid)%lon1)
  istat = nf90_put_att(ncid, NF90_GLOBAL, 'lat0', idx_nc(mytid)%lat0)
  istat = nf90_put_att(ncid, NF90_GLOBAL, 'lat1', idx_nc(mytid)%lat1)
  ! the raw %data time level indices, so the mapping stays recoverable
  istat = nf90_put_att(ncid, NF90_GLOBAL, 'itp', itp)
  istat = nf90_put_att(ncid, NF90_GLOBAL, 'itc', itc)
  call add_creation_time(ncid)

  ! coordinates are indexed in the glon/glat space -> idx_nc
  call init_spacial_dims_subdomain(ncid, &
         idx_nc(mytid)%lon0, idx_nc(mytid)%lon1, &
         idx_nc(mytid)%lat0, idx_nc(mytid)%lat1, &
         dim_ids, coord_var_ids)

  ! define one variable per field ------------------------------------

  allocate(var_ids(size(f4d)))

  do i = 1, size(f4d)
    ! dim_ids is (lon, lat, lev, ilev, timelevel)
    select case(trim(f4d(i)%vcoord))
      case('interfaces')
        lev_dim = dim_ids(4)
      case('midpoints')
        lev_dim = dim_ids(3)
      case default
        lev_dim = dim_ids(3)
        write(*,*) 'panic_write: unknown vcoord <',trim(f4d(i)%vcoord), &
                   '> for ',trim(f4d(i)%short_name),' - assuming midpoints'
    end select

    istat = nf90_def_var(ncid, trim(f4d(i)%short_name), NF90_DOUBLE, &
                         (/dim_ids(1), dim_ids(2), lev_dim, dim_ids(5)/), var_ids(i))
    if (istat /= NF90_NOERR) call handle_ncerr(istat, &
      'panic_write: Error defining variable '//trim(f4d(i)%short_name))

    istat = nf90_put_att(ncid, var_ids(i), "long_name", trim(f4d(i)%long_name))
    units = f4d(i)%units
    if(len_trim(units)==0) units = ' '   ! never write a zero length attribute
    istat = nf90_put_att(ncid, var_ids(i), "units", trim(units))
    istat = nf90_put_att(ncid, var_ids(i), "vcoord", trim(f4d(i)%vcoord))
  end do

  ! also include diagnostig 3d fields
  allocate(var_ids_3d(size(f3d)))

  do i = 1, size(f3d)
    select case(trim(f3d(i)%vcoord))
      case('interfaces')
        lev_dim = dim_ids(4)
      case('midpoints')
        lev_dim = dim_ids(3)
      case default
        lev_dim = dim_ids(3)
        write(*,*) 'panic_write: unknown vcoord <',trim(f3d(i)%vcoord), &
                   '> for ',trim(f3d(i)%short_name),' - assuming midpoints'
    end select

    istat = nf90_def_var(ncid, trim(f3d(i)%short_name), NF90_DOUBLE, &
                         (/dim_ids(1), dim_ids(2), lev_dim/), var_ids_3d(i))
    if (istat /= NF90_NOERR) call handle_ncerr(istat, &
      'panic_write: Error defining variable '//trim(f3d(i)%short_name))

    istat = nf90_put_att(ncid, var_ids_3d(i), "long_name", trim(f3d(i)%long_name))
    units = f3d(i)%units
    if(len_trim(units)==0) units = ' '   ! never write a zero length attribute
    istat = nf90_put_att(ncid, var_ids_3d(i), "units", trim(units))
    istat = nf90_put_att(ncid, var_ids_3d(i), "vcoord", trim(f3d(i)%vcoord))
  end do

  istat = nf90_enddef(ncid)
  if (istat /= NF90_NOERR) call handle_ncerr(istat,'panic_write: Error leaving define mode')

  ! write values -----------------------------------------------------

  call write_spacial_dims_subdomain(ncid, &
         idx_nc(mytid)%lon0, idx_nc(mytid)%lon1, &
         idx_nc(mytid)%lat0, idx_nc(mytid)%lat1, &
         coord_var_ids)

  allocate(field_reshaped(nlons,nlats,nlevp1), stat=istat)
  if(istat/=0) then
    write(*,*) 'panic_write: could not allocate work array - closing file'
    istat = nf90_close(ncid)
    return
  end if

  count_p = (/nlons, nlats, nlevp1, 1/)
  start_p(1:3) = 1

  do i = 1, size(f4d)
    do k = 1, 2
      ! internal storage is (lev,lon,lat), the file convention is (lon,lat,lev),
      ! cf. nc_tiegcm_write in result_file_writer.F90.
      ! %data is sliced with idx_intern (TIE-GCM internal index space)
      field_reshaped = reshape( &
        f4d(i)%data(:, &
                    idx_intern(mytid)%lon0:idx_intern(mytid)%lon1, &
                    idx_intern(mytid)%lat0:idx_intern(mytid)%lat1, &
                    time_slot(k)), &
        shape=(/nlons,nlats,nlevp1/), order=(/3,1,2/))

      start_p(4) = k
      istat = nf90_put_var(ncid, var_ids(i), field_reshaped, &
                           start=start_p, count=count_p)
      if (istat /= NF90_NOERR) call handle_ncerr(istat, &
        'panic_write: Error writing '//trim(f4d(i)%short_name))
    end do
  end do

  do i = 1, size(f3d)
    field_reshaped = reshape( &
      f3d(i)%data(:, &
                  idx_intern(mytid)%lon0:idx_intern(mytid)%lon1, &
                  idx_intern(mytid)%lat0:idx_intern(mytid)%lat1), &
      shape=(/nlons,nlats,nlevp1/), order=(/3,1,2/))

    istat = nf90_put_var(ncid, var_ids_3d(i), field_reshaped)
    if (istat /= NF90_NOERR) call handle_ncerr(istat, &
      'panic_write: Error writing '//trim(f3d(i)%short_name))
  end do

  deallocate(field_reshaped)
  deallocate(var_ids)
  deallocate(var_ids_3d)

  istat = nf90_close(ncid)
  if (istat /= NF90_NOERR) call handle_ncerr(istat,'panic_write: Error closing '//trim(filename))

  write(*,*) 'panic_write: wrote ', trim(filename)

end subroutine panic_write

!> Creates the instant-debug NetCDF file (on rank 0) with ensemble/time/spatial dimensions and global attributes.
subroutine init_nc_debug(name)

    USE mod_assimilation, only: dim_ens   

    use netcdf_functionality, only: init_spacial_dims, init_temporal_dims, add_creation_time

    implicit none

    character(len=*), intent(in) :: name

    ! local
    character(len=80) :: char80

    integer :: dim_id_ensemble

    ncfile = trim(name)

    if ( rank_world == 0) then

        istat = nf_create(ncfile, NF_NETCDF4, ncid)
        if (istat /= NF_NOERR) call handle_ncerr(istat, 'Error creating '//ncfile)

        ! write global attributes
        char80 = 'TIE-GCM debug output'

        istat = NF_PUT_ATT_TEXT(ncid, NF_GLOBAL, 'title', LEN_TRIM(char80),  TRIM(char80))
        if (istat /= NF_NOERR) call handle_ncerr(istat,'Error writing global attributte')

        call add_creation_time(ncid)

        ! Define ensemble dimensions
        istat = NF_DEF_DIM(ncid, 'ensemble',  dim_ens, dim_id_ensemble)
        if (istat /= NF_NOERR) call handle_ncerr(istat,'Error defining dim step')

        ! add temporal dimension
        call init_temporal_dims(ncid)

        ! add spatial dimensions
        call init_spacial_dims(ncid)

        istat = nf_close(ncid)
        if (istat /= NF_NOERR) call handle_ncerr(istat,'Error closing '//ncfile)

    end if

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

end subroutine init_nc_debug

!> Advances the instant-debug NetCDF file's unlimited time index and writes the current model time.
subroutine advance_nc_debug

    use netcdf_functionality, only: add_model_time

    implicit none

    ulimlen = ulimlen + 1

    if ( rank_world == 0) then

        istat = nf_open(ncfile, NF_WRITE, ncid)

        call add_model_time(ncid, ulimlen, .false.)

        istat = nf_close(ncid)
        if (istat /= NF_NOERR) call handle_ncerr(istat,'Error closing '//ncfile)

    end if

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

end subroutine advance_nc_debug

!> Immediately writes one 3D field to the instant-debug NetCDF file for the current ensemble member and time step, defining the variable on first use.
!!
!! vcoord is either 'mid' or 'int'.
subroutine instant_write( var_name, vcoord, f3d)

    ! tie-gcm
    use params_module, only: nlevp1
    use fields_module, only:  levd0,levd1, lond0,lond1, latd0,latd1
    use mpi_module, only: mytid

    ! intern
    use state_module, only: idx_nc, idx_intern
    use mod_parallel_pdaf, only: task_id

    implicit none

    !!! important lesson learned here
    !  explicitly specify bounds! 
    !
    !  real, dimension(levd0:levd1, lond0:lond1, latd0:latd1, 2) :: test
    !  call fun( test(:,:,:,1) )
    !  explicitly specified bounds are lost when entering function !
    !
    ! Thus, we need to define
    ! 
    ! real, dimension(levd0:levd1,lond0:lond1,latd0:latd1), intent(in) :: f3d
    !
    ! real, dimension(:,:,:), intent(in) :: will have indices staring at 1
    !

    real, dimension(levd0:levd1,lond0:lond1,latd0:latd1), intent(in) :: f3d
    character(len=*), intent(in) :: var_name
    character(len=3), intent(in) :: vcoord

    integer :: var_idx
    integer :: lev_dim

    integer, dimension(5) :: start_p, count_p

    integer :: dim_id_ensemble, dim_id_time, &
               dim_id_lon, dim_id_lat, dim_id_lev, dim_id_ilev

    istat = nf_open_par(ncfile, NF_WRITE, MPI_COMM_WORLD,  MPI_INFO_NULL, ncid)
    if (istat /= NF_NOERR) call handle_ncerr(istat,'Error opening parallel netcdf '// ncfile)

    istat = nf_inq_unlimdim(ncid, dim_id_time)
    if (istat /= NF_NOERR) call handle_ncerr(istat,'Error INQUIRE unlimited dim')

    istat = nf_inq_dimid(ncid, "ensemble", dim_id_ensemble)
    istat = nf_inq_dimid(ncid, "lon", dim_id_lon)
    istat = nf_inq_dimid(ncid, "lat", dim_id_lat)
    istat = nf_inq_dimid(ncid, "lev", dim_id_lev)
    istat = nf_inq_dimid(ncid, "ilev", dim_id_ilev)

    istat = nf_inq_varid(ncid, "mtime", var_id_mtime)
    if (istat /= NF_NOERR) call handle_ncerr(istat, 'Error INQUIRE mtime')

    istat = nf_inq_varid(ncid, "time", var_id_time)
    if (istat /= NF_NOERR) call handle_ncerr(istat, 'Error INQUIRE time')

    if (var_list_len < 1) then
        var_idx = 0
    else
        var_idx = findloc( var_list(1:var_list_len), var_name, dim=1)
    end if

    ! create new varibale if not existing
    if(  var_idx < 1 ) then
        var_list_len = var_list_len + 1

        var_idx = var_list_len

        if( var_list_len > maxDebugFields ) then
            call shutdown( 'to many debug fields. increase maxDebugFields in debugFun.F ' )
        end if

        write(var_list(var_idx), '(a)' ) trim(var_name)
        istat = nf_redef(ncid)
        if (istat /= NF_NOERR) call handle_ncerr(istat,'Error entering define mode')

        select case (trim(vcoord) )
            case ('int')
                lev_dim = dim_id_ilev
            case ('mid')
                lev_dim = dim_id_lev
            case default
                lev_dim = -1
        end select

        istat = NF_DEF_VAR(ncid, var_list(var_idx), NF_DOUBLE, 5, &
                                (/dim_id_ensemble, lev_dim, dim_id_lon, dim_id_lat, dim_id_time/), &
                                var_id(var_idx) )
        if (istat /= NF_NOERR) call handle_ncerr(istat,'Error defining variable '// var_list(var_idx) )

        istat = NF_ENDDEF(ncid) 
        if (istat /= NF_NOERR) call handle_ncerr(istat,'Error leaving define mode')
    else

        istat = nf_inq_varid (ncid, var_list(var_idx), var_id(var_idx));
        if (istat /= NF_NOERR) call handle_ncerr(istat, 'Error inq '//var_list(var_idx))
    end if

    ! write data
    istat = nf_var_par_access(ncid, var_id(var_idx), NF_COLLECTIVE);
    if (istat /= NF_NOERR) call handle_ncerr(istat, 'Error setting par_access'//var_list(var_idx))

    start_p(1) = task_id ! ensemble
    start_p(2) = 1       ! level
    start_p(3) = idx_nc(mytid)%lon0
    start_p(4) = idx_nc(mytid)%lat0
    start_p(5) = ulimlen

    count_p(1) = 1
    count_p(2) = nlevp1
    count_p(3) = idx_nc(mytid)%nlons
    count_p(4) = idx_nc(mytid)%nlats
    count_p(5) = 1

    istat = nf_put_vara_double(ncid, var_id(var_idx), start_p, count_p,  &
                f3d(                      :,                             &
                    idx_intern(mytid)%lon0:idx_intern(mytid)%lon1,       &
                    idx_intern(mytid)%lat0:idx_intern(mytid)%lat1))
    if (istat /= NF_NOERR) call handle_ncerr(istat,'Error putting '//var_name)

    istat = nf_close(ncid)

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
end subroutine instant_write

end module nc_debug
