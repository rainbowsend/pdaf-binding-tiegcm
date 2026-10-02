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
! Defines regular-grid observation datasets (ESMF grid, NetCDF I/O) and maps them onto TIE-GCM grid cells/state vector.

module grid_observation_module
  ! observation located on a grid with potentially several parameters per grid cell

  ! extern
  use esmf
  use mpi_f08

  ! intern
  use mapping_module, only: mapping

  implicit none

  type :: field_bundle
     type( ESMF_Field ), dimension(:), allocatable :: bundle
     character(len=16), dimension(:), allocatable :: names
     contains
     procedure, pass(this), public :: field => field_bundle_get_field
     procedure, pass(this), public :: size => field_bundle_get_number_fields
     procedure, pass(this), public :: destroy => field_bundle_destroy
  end type field_bundle

  ! Two classes are used to represent netcdf file holding observations on a regular grid (output of tme::evaluate)
  !
  ! reg_grid_dataset_root holds the data that is common to all groups,i.e., spatial and temporal inforamtion
  ! reg_grid_dataset_group represents a group of the netcdf file with fields build on the dimensions of the root
  !

  type reg_grid_dataset_root
    integer :: ncid
    character(len=128) :: nc_file
    ! full coordinate vector
    real, allocatable, dimension(:) :: lon, lat, alt
    type(ESMF_Time), allocatable, dimension(:) :: epochs
    integer :: alt_first, alt_last
    ! spatial grid related variables (subdomain)
    integer :: nlon       ! lon subdomain
    integer :: nlat       ! lat subdomain
    integer :: nalt       ! alt subdomain
    integer :: n          ! nlon * nlat * nalt
    integer :: offset_lon ! offsets at this rank
    integer :: offset_lat !
    type( ESMF_Grid ) :: grid

    integer :: t_idx = 1  ! current time index

    character(len=80) :: first_timestr, last_timestr

    contains
    procedure, pass(this) :: destroy => reg_grid_dataset_root_destroy
    procedure, pass(this) :: add_spatial_dim =>reg_grid_dataset_root_add_spatial_dim
  end type reg_grid_dataset_root

  type reg_grid_dataset_group
    type(reg_grid_dataset_root), pointer :: root => NULL()
    integer :: grp_ncid
    type(mapping) map
    ! field related variabels
    type(field_bundle) :: bundle
    !
    contains
    procedure, pass(this) :: destroy => reg_grid_dataset_group_destroy
    procedure, pass(this) :: read_current_epoch => reg_grid_dataset_group_read_current_epoch
  end type reg_grid_dataset_group

  contains

  !> Creates a field bundle with one new ESMF_Field per given name, all defined on the given ESMF grid.
  function construct_field_bundle(grid, names) result(this)

    implicit none

    ! arguments
    type( ESMF_Grid ), intent(in) :: grid
    character(len=*), dimension(:), intent(in) :: names

    ! returns
    type (field_bundle) :: this

    ! local
    integer :: i, rc

    ! ATTENTION allocate(,source=) not working here. Leads to run time error
    ! probably by reason of character(len=*)
    allocate(this%names(size(names)))
    do i=1,size(names)
        write(this%names(i),'(a)') names(i)
    end do

    allocate(this%bundle(size(names)))

    do i=1,size(names)
      this%bundle(i) = ESMF_FieldCreate( grid, &
                                         typekind = ESMF_TYPEKIND_R8, &
                                         name = trim(names(i)),&
                                         rc=rc)
      if(ESMF_LogFoundError(rc,msg="construct_field_bundle:ESMF_FieldCreate", &
        rcToReturn=rc)) then
        call shutdown('ESMF error in construct_field_bundle:ESMF_FieldCreate')
      end if
    end do

  end function

  !> Returns the number of fields in the bundle.
  function field_bundle_get_number_fields(this) result(n)

    implicit none

    ! arguments
    class (field_bundle), intent(inout) :: this

    ! returns
    integer :: n

    n = size(this%bundle)
  end function


  !> Returns a pointer to a field in the bundle, looked up by name or by index.
  function field_bundle_get_field(this,name,idx) result(field)

    implicit none

    ! arguments
    class (field_bundle), intent(inout) :: this
    character(len=*), intent(in), optional :: name
    integer, intent(in), optional :: idx

    ! result
    type( ESMF_Field ), pointer :: field

    ! local
    integer :: idx_

    idx_=0

    if(present(name)) then
      idx_ = findloc( this%names, name, dim=1 )
      if ( idx_ < 1) then
        call shutdown('field_bundle_get_field: field '// name // ' not found in field bundle')
      end if
    else if(present(idx)) then
      if((idx > size(this%names)) .or. (idx < 1)) then
        call shutdown('field_bundle_get_field: idx is out of bounds')
      end if
      idx_=idx
    else
      call shutdown('field_bundle_get_field: invalid call of field_bundle_get_field')
    end if

    if ( idx_ > 0) then
      field => field_ptr(this%bundle(idx_))
    else
      field => null()
    end if

  end function

  !> Destroys all ESMF fields in the bundle and deallocates its arrays.
  subroutine field_bundle_destroy(this)

    implicit none

    ! arguments
    class(field_bundle), intent(inout) :: this

    ! local
    integer :: i, rc

    if (allocated(this%bundle)) then
      do i = 1,size(this%bundle)
         if(ESMF_FieldIsCreated(this%bundle(i),rc=rc))then
            call ESMF_FieldDestroy( this%bundle(i), noGarbage=.true.,rc=rc)
         end if
      end do
      deallocate(this%bundle)
    end if

    if (allocated(this%names)) deallocate(this%names)

  end subroutine

  !> Returns a pointer targeting the given ESMF_Field.
  function field_ptr(field) result(ptr)
    type( ESMF_Field ), target, intent(in) :: field
    type( ESMF_Field ), pointer :: ptr
    ptr=>field
  end function

  !> Opens a regular-grid observation NetCDF file, reads its coordinates and epochs, and builds the shared ESMF grid for the dataset.
  function construct_reg_grid_dataset_root(nc_file,altmin, altmax) result (this)

    ! extern
    use netcdf

    ! tie-gcm
    use mpi_module, only: TIEGCM_WORLD
    use nchist_module, only: handle_ncerr

    implicit none

    ! arguments
    character(len=*), intent(in) :: nc_file
    real, intent(in), optional :: altmin, altmax

    ! returns
    type(reg_grid_dataset_root) :: this

    ! local
    logical :: there
    integer :: istat

    inquire( file=nc_file, exist=there )
    if ( there .neqv. .true. ) then
        call shutdown('construct_reg_grid_dataset_group: netcdf file ' // nc_file // ' does not exist')
    end if

    write(this%nc_file,'(a)') trim(nc_file)

    istat = nf90_open(path=nc_file, &
                      mode=NF90_NOWRITE, &
                      ncid=this%ncid, &
                      comm= TIEGCM_WORLD%MPI_VAL, &
                      info=MPI_INFO_NULL%MPI_VAL)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error opening netcdf' // this%nc_file)

    ! reads coordinate vars (grid dimensions)
    call read_reg_grid_nc( this%ncid, this%alt, this%lon, this%lat, this%epochs )

    call ESMF_TimeGet( this%epochs(1),&
                       timeString=this%first_timestr)
    call ESMF_TimeGet( this%epochs(ubound(this%epochs,1)),&
                       timeString=this%last_timestr)

    if(present(altmin))then
      this%alt_first = minloc(abs(this%alt-altmin),dim=1,back=.true.)
    else
      this%alt_first = 1
    end if

    if(present(altmax))then
      this%alt_last = minloc(abs(this%alt-altmax),dim=1)
    else
      this%alt_last = size(this%alt)
    end if

    ! create esmf grid with given coordinates
    call reg_grid_dataset_group_init_ESMF_grid( this,&
                          this%alt(this%alt_first:this%alt_last),&
                          this%lon,&
                          this%lat)


  end function construct_reg_grid_dataset_root

  !> Closes the dataset's NetCDF file and deallocates its coordinate arrays and ESMF grid.
  subroutine reg_grid_dataset_root_destroy(this)

    ! extern
    use netcdf

    ! tie-gcm
    use nchist_module, only: handle_ncerr

    implicit none

    ! arguments
    class(reg_grid_dataset_root) :: this

    ! local
    integer :: istat

    ! close netcdf files
    istat = nf90_close(this%ncid)

    if ( allocated(this%epochs)) deallocate(this%epochs)
    if ( allocated(this%lon)) deallocate(this%lon)
    if ( allocated(this%lat)) deallocate(this%lat)
    if ( allocated(this%alt)) deallocate(this%alt)

    if(ESMF_GridIsCreated(this%grid))then
      call ESMF_GridDestroy(this%grid, noGarbage=.true.)
    end if

  end subroutine reg_grid_dataset_root_destroy

  !> Writes this dataset's spatial (lon/lat/alt) dimensions and coordinate variables into another already-open NetCDF file.
  subroutine reg_grid_dataset_root_add_spatial_dim(this,ncid)

    ! extern
    use netcdf

    ! tie-gcm
    use nchist_module, only: handle_ncerr

    implicit none

    ! arguments
    class(reg_grid_dataset_root), intent(in) :: this
    integer, intent(in):: ncid

    ! local
    integer :: istat
    integer :: dim_id_alt, dim_id_lat, dim_id_lon
    integer :: var_id_alt, var_id_lon, var_id_lat

    istat = nf90_def_dim(ncid, "alt", size(this%alt), dim_id_alt)
    istat = nf90_def_dim(ncid, "lat", size(this%lat), dim_id_lat)
    istat = nf90_def_dim(ncid, "lon", size(this%lon), dim_id_lon)

    istat = nf90_def_var(ncid, "alt", NF90_DOUBLE, dim_id_alt, var_id_alt )
    istat = nf90_def_var(ncid, "lat", NF90_DOUBLE, dim_id_lat, var_id_lat )
    istat = nf90_def_var(ncid, "lon", NF90_DOUBLE, dim_id_lon, var_id_lon )

    istat = nf90_put_att(ncid, var_id_alt, "units", "m")
    istat = nf90_put_att(ncid, var_id_lat, "units", "degrees_north")
    istat = nf90_put_att(ncid, var_id_lon, "units", "degrees_east")

    istat = nf90_put_var(ncid, var_id_alt, this%alt )
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error putting alt ')
    istat = nf90_put_var(ncid, var_id_lat, this%lat )
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error putting lat ')
    istat = nf90_put_var(ncid, var_id_lon, this%lon )
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error putting lon ')

  end subroutine

  !> Opens a group in the observation NetCDF file, creates its field bundle, and builds the process-to-grid index mapping.
  function construct_reg_grid_dataset_group(root, grp_name, field_names, grid_name) result (this)

    ! extern
    use netcdf

    ! tie-gcm
    use mpi_module, only: TIEGCM_WORLD, mytid, ntask
    use nchist_module, only: handle_ncerr

    ! intern
    use mapping_module, only: construct_mapping, print_mapping

    implicit none

    ! arguments
    type(reg_grid_dataset_root), intent(in), target :: root
    character(len=*), intent(in), optional :: grp_name
    character(len=*), dimension(:), intent(in) :: field_names
    character(len=*), intent(in) :: grid_name

    ! returns
    type(reg_grid_dataset_group) :: this

    ! local
    integer :: i
    integer :: istat

    integer, dimension(:,:), allocatable :: F_R_size
    integer :: ierr
    integer, dimension(ntask) :: subdomain_size

    this%root => root

    if (present(grp_name)) then
      istat = nf90_inq_ncid( this%root%ncid, trim(grp_name), this%grp_ncid )
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error getting group id for group '// grp_name)
    else
      this%grp_ncid = this%root%ncid
    end if

    this%bundle = construct_field_bundle(this%root%grid, field_names)

    ! compute mapping
    call MPI_Allgather( this%root%n, 1, MPI_Integer, &
                        subdomain_size, 1, MPI_Integer, &
                        TIEGCM_WORLD, ierr)

    allocate(F_R_size(1:this%bundle%size(), 0:ntask-1) )

    ! All fields have the same size
    do i = 0, ntask-1
      F_R_size(:,i) = subdomain_size(i+1)
    end do

    this%map = construct_mapping(F_R_size, field_names, grid_name)
    deallocate(F_R_size)
    if(mytid==0) call print_mapping(this%map)


  end function construct_reg_grid_dataset_group

  !> Destroys the group's field bundle and its process-to-grid index mapping.
  subroutine reg_grid_dataset_group_destroy( this )

    ! intern
    use mapping_module, only: deconstruct_mapping

    implicit none

    ! arguments
    class(reg_grid_dataset_group), intent(inout) :: this

    call this%bundle%destroy
    call deconstruct_mapping(this%map)

  end subroutine reg_grid_dataset_group_destroy

  !> Distributes the regular grid across processes and creates the shared ESMF grid with per-process coordinate subdomains.
  subroutine reg_grid_dataset_group_init_ESMF_grid(this, lev, lon, lat )

    ! tie-gcm
    use mpi_module, only: TIEGCM_WORLD, ntaski, ntaskj, mytid

    ! intern
    use array_print_module, only: printMat
    use tiegcm_optimized_interpolator, only: distriubute_grid

    ! arguments
    type (reg_grid_dataset_root),  intent(inout) :: this
    real, dimension(:), intent(in)  :: lev, lon, lat ! full (gathered) coordinates

    ! local
    integer :: rc
    real(ESMF_KIND_R8), pointer :: coordLon(:), coordLat(:), coordLev(:)
    integer :: lbnd(1),ubnd(1)
    integer :: nlons_task(ntaski,ntaskj)
    integer :: nlats_task(ntaski,ntaskj)
    integer :: ierr

    integer :: count_lon, count_lat

    ! TODO check this
    call distriubute_grid(lon, lat, this%offset_lon, count_lon, this%offset_lat, count_lat)

    call MPI_Allgather( &
      count_lon,  1, MPI_Integer, &
      nlons_task, 1, MPI_Integer, &
    TIEGCM_WORLD, ierr)

!     call printMat(nlons_task,'lons')

    call MPI_Allgather( &
      count_lat,  1, MPI_Integer, &
      nlats_task, 1, MPI_Integer, &
    TIEGCM_WORLD, ierr)

!     call printMat(nlats_task,'lats')


    this%nlon = count_lon
    this%nlat = count_lat
    this%nalt = size(lev)
    this%n = this%nlon * this%nlat * this%nalt

    write(*,'(a)') 'read grid information: Dimension of Grid'
    write(*,'(a,i4,a,/,a,i4,/,a,i4,/,a,i4,/,a,i8)') &
               'subdomain (model rank:', mytid, ')', &
               '   lon: ', this%nlon, &
               '   lat: ', this%nlat, &
               '   alt: ', this%nalt, &
               ' total: ', this%n
    write(*,'(a,/,a,i4,/,a,i4,/,a,i4,/,a,i8)') &
               'full grid', &
               '   lon: ', size(lon), &
               '   lat: ', size(lat), &
               '   alt: ', size(lev), &
               ' total: ', size(lon) * size(lat) * size(lev)

    this%grid = ESMF_GridCreate1PeriDim( &
         coordSys = ESMF_COORDSYS_SPH_DEG, &
         countsPerDEDim1=nlons_task(:,1), coordDep1=(/1/), &
         countsPerDEDim2=nlats_task(1,:), coordDep2=(/2/), &
         countsPerDEDim3=(/size(lev)/), coordDep3=(/3/), &
         indexflag=ESMF_INDEX_GLOBAL, &
         minIndex=(/1,1,1/),rc=rc)

    if(ESMF_LogFoundError(rc,msg="read_reg_grid_nc:ESMF_GridCreate1PeriDim", &
        rcToReturn=rc)) then
        call shutdown('ESMF error in read_reg_grid_nc:ESMF_GridCreate1PeriDim')
    end if

    call ESMF_GridAddCoord( this%grid, &
         staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)

    if(ESMF_LogFoundError(rc,msg="read_reg_grid_nc:ESMF_GridAddCoord", &
        rcToReturn=rc)) then
        call shutdown('ESMF error in read_reg_grid_nc:ESMF_GridAddCoord')
    end if

    call ESMF_GridGetCoord( this%grid, coordDim=1, localDE=0, &
       computationalLBound=lbnd, computationalUBound=ubnd, &
       farrayPtr=coordLon, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)

    if(ESMF_LogFoundError(rc,msg="read_reg_grid_nc:ESMF_GridGetCoord", &
       rcToReturn=rc)) then
        call shutdown('ESMF error in read_reg_grid_nc:ESMF_GridGetCoord')
    end if

    call ESMF_GridGetCoord( this%grid, coordDim=2, localDE=0, &
       computationalLBound=lbnd, computationalUBound=ubnd, &
       farrayPtr=coordLat, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)

    if(ESMF_LogFoundError(rc,msg="read_reg_grid_nc:ESMF_GridGetCoord", &
       rcToReturn=rc)) then
        call shutdown('ESMF error in read_reg_grid_nc:ESMF_GridGetCoord')
    end if

    call ESMF_GridGetCoord( this%grid, coordDim=3, localDE=0, &
       computationalLBound=lbnd, computationalUBound=ubnd, &
       farrayPtr=coordLev, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)

    if(ESMF_LogFoundError(rc,msg="read_reg_grid_nc:ESMF_GridGetCoord", &
       rcToReturn=rc)) then
        call shutdown('ESMF error in read_reg_grid_nc:ESMF_GridGetCoord')
    end if

    coordLon( this%offset_lon : this%offset_lon + count_lon -1 ) = &
         lon( this%offset_lon : this%offset_lon + count_lon -1 )

    coordLat( this%offset_lat : this%offset_lat + count_lat -1 ) = &
         lat( this%offset_lat : this%offset_lat + count_lat -1 )

    coordLev = lev
  end subroutine reg_grid_dataset_group_init_ESMF_grid

  !> Reads this group's fields at the current time index (t_idx) from the NetCDF file into the ESMF field bundle.
  subroutine reg_grid_dataset_group_read_current_epoch( this )

    ! extern
    use netcdf

    ! tie-gcm
    use nchist_module, only: handle_ncerr

    ! intern
    use netcdf_functionality, only : orderNetcdfVar

    implicit none

    ! arguments
    class (reg_grid_dataset_group),  intent(inout) :: this

    ! local
    integer :: var_id

    integer :: startp(4)
    integer :: countp(4)
    integer :: rc, istat

    integer :: field_id, i

    integer :: startp2(4)
    integer :: countp2(4)

    integer, dimension(4)  :: order

    real, allocatable, dimension(:,:,:) :: tmp_data

    real(ESMF_KIND_R8), pointer :: field_data(:,:,:)

    character(len=16) :: field_name
    character(len=512) :: buff

!     integer, save :: counter = 1

    startp = 1
    countp = 1
    startp2 = 1
    countp2 = 1

!     write(*,'(a,i6)') 'loading the fields: into esmf field.', this%map%vector_name, '. Time index of observation in netcdf file: ',  this%root%t_idx

     write(buff,*) this%bundle%names
     write(*,*) 'loading the fields:', trim(buff), ' at time index ',  this%root%t_idx,&
                ' from netcdf file: ', trim(this%root%nc_file), ' into ESMF Field'

    startp(4) = this%root%t_idx

    call ESMF_GridGetFieldBounds(this%root%grid, &
        localDe     = 0          , &
        totalCount  = countp(1:3), &
        totalLBound = startp(1:3), &
        rc          = rc             )

    startp(3)=this%root%alt_first

!!! DEBUG
!     write(*,'(a,I4,/,a,4I10,/,a,4I10,/,a,4I10)') &
!                 'rank ', mytid,  &
!                 'start   ', startp, &
!                 'count   ', countp, &
!                 'end     ', startp+countp-1
!!! DEBUG

    fieldLoop: do field_id = 1, this%bundle%size()

        call ESMF_FieldGet( this%bundle%bundle(field_id), name=field_name, rc=rc )

        istat = nf90_inq_varid( this%grp_ncid, trim(field_name), var_id)
        if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error getting var id ' // field_name)

        ! We need the dimensions in order 'lon ', 'lat ', 'lev ', 'time'
        ! the following call will return the order in the netcdf file
        call orderNetcdfVar( this%grp_ncid, var_id, (/'lon ', 'lat ', 'alt ', 'time'/), order)

        ! startp and countp assume the order specified above
        ! We need to rearange them to be compatible to the netcdf file
        do i = 1,4
            startp2(i) = startp( order(i) )
            countp2(i) = countp( order(i) )
        end do

    !!! DEBUG
    !     write(*,'(a,I4,/,a,4I10,/,a,4I10,/,a,4I10)') &
    !             'rank ', mytid,  &
    !             'start   ', startp2, &
    !             'count   ', countp2, &
    !             'end     ', startp2+countp2-1
    !!! DEBUG

        allocate (tmp_data( countp(order(1)), countp(order(2)) , countp(order(3))   ))

        ! read netcdf directly into ESMF field is possibel but we need to rearange netcdf first
        istat = nf90_get_var(ncid=this%grp_ncid,&
                             varid=var_id,&
                             values=tmp_data,&
                             start=startp2,&
                             count=countp2)
        if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error getting variable ' // field_name)

        ! Get the pointer of the field
        call ESMF_FieldGet( this%bundle%bundle(field_id), localDe=0, farrayPtr=field_data, &
            rc=rc)

        ! reshape netcdf data to correct order
        field_data = reshape( tmp_data, countp(1:3), order=order(1:3) )

!         slicedim = 3
!         call printMat(field_data, field_name)

        call ESMF_FieldValidate( this%bundle%bundle(field_id), rc=rc)
        if(ESMF_LogFoundError(rc, msg="reg_grid_dataset_group_read_current_epoch:ESMF_FieldValidate", &
            rcToReturn=rc)) then
            call shutdown("reg_grid_dataset_group_read_current_epoch:ESMF_FieldValidate")
        end if

        deallocate( tmp_data )

! ! write field to check wether reading routine is correct
!         call ESMF_FieldWrite( this%fields(field_id), &
!             fileName="grid_obs_in.nc",   &
!             variableName=field_name, &
!             iofmt=ESMF_IOFMT_NETCDF, &
!             timeslice = counter, &
!             overwrite = .true., &
!             rc=rc)

    end do fieldLoop
!     counter =counter + 1
  end subroutine reg_grid_dataset_group_read_current_epoch

  !> Reads the coordinate (lev/lon/lat) and time variables of a regular-grid observation NetCDF file.
  subroutine read_reg_grid_nc(ncid, lev, lon, lat, time)

    ! extern
    use netcdf

    ! tie-gcm
    use nchist_module, only: handle_ncerr

    ! intern
    use netcdf_functionality, only: read_coordinate_var
    use time_module, only: start_epoch_form_att

    implicit none

    ! arguments
    integer, intent(in) :: ncid
    real, dimension(:), allocatable, intent(inout) :: lev, lon, lat
    type(ESMF_Time), allocatable, dimension(:), intent(inout) :: time

    ! local
    real, dimension(:), allocatable :: real_time

    integer :: istat
    integer :: rc
    integer :: i

    integer :: var_id_lev, var_id_lon, var_id_lat, var_id_time

    character(len=80) :: time_unit
    type(ESMF_Time) :: start_epoch
    type(ESMF_TimeInterval) :: timeinterval


    call read_coordinate_var(ncid,'alt',var_id_lev,lev)
    call read_coordinate_var(ncid,'lon',var_id_lon,lon)
    call read_coordinate_var(ncid,'lat',var_id_lat,lat)
    call read_coordinate_var(ncid,'time',var_id_time,real_time)

    ALLOCATE(time(size(real_time)))

    istat =nf90_get_att(ncid, var_id_time, 'units', time_unit)

    call start_epoch_form_att( trim(time_unit), start_epoch )

    do i = 1, size(real_time)

        call ESMF_TimeIntervalSet(timeinterval, s_r8=real_time(i), rc=rc)
        if(ESMF_LogFoundError(rc,msg="read_reg_grid_nc:ESMF_TimeIntervalSet", &
            rcToReturn=rc)) then
            call shutdown("read_reg_grid_nc:ESMF_TimeIntervalSet")
        end if

        time(i) = start_epoch + timeinterval
!!! DEBUG time parsing
!         if(mytid==0) call ESMF_TimePrint(this%epochs(i), options="string", rc=rc)
    end do

    deallocate(real_time)

  end subroutine read_reg_grid_nc

  !> Fills the given array with synthetic test values evaluated at the observation grid's coordinates.
  subroutine fill_synthetic_obs_p( obsgrid, out )

    ! intern
    use test_values, only: synthetic_grid

    implicit none

    ! arg
    type(ESMF_grid), intent(in) :: obsgrid
    real(ESMF_KIND_R8), dimension(:,:,:), intent(out) :: out

    ! local
    integer :: rc
    logical :: LogFoundError
    real(ESMF_KIND_R8), pointer :: coordLon(:), coordLat(:), coordLev(:)

    ! Get pointer to coordinate coords
    call ESMF_GridGetCoord(obsgrid, coordDim=1, localDE=0, &
        farrayPtr=coordLon, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    LogFoundError = ESMF_LogFoundError(rc, msg="fill_synthetic_obs_p ESMF_GridGetCoord 1", rcToReturn=rc)

    call ESMF_GridGetCoord(obsgrid, coordDim=2, localDE=0, &
        farrayPtr=coordLat, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    LogFoundError = ESMF_LogFoundError(rc, msg="fill_synthetic_obs_p ESMF_GridGetCoord 2", rcToReturn=rc)

    call ESMF_GridGetCoord(obsgrid, coordDim=3, localDE=0, &
        farrayPtr=coordLev, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    LogFoundError = ESMF_LogFoundError(rc, msg="fill_synthetic_obs_p ESMF_GridGetCoord 3", rcToReturn=rc)

    out = synthetic_grid(coordLon, coordLat, coordLev)

  end subroutine fill_synthetic_obs_p

  !> Fills the observation coordinate array with each grid point's (lon, lat, alt), in flattened field order.
  subroutine get_ocoord_lon_lat_lev(ocoord_p,obs)
    ! ocoord_p(1) longitude in radian
    ! ocoord_p(2) latitude  in radian
    ! ocoord_p(3) altitude  in meters, only if requested (see below)
    !
    ! ocoord_p is allocated with one row per coordinate used for the distance
    ! computation (PDAF's ncoord). Quantities without vertical extent, i.e.
    ! VTEC, are localized horizontally only and therefore provide two rows.
    ! The altitude is written only if there is a row for it.

    ! extern
    use esmf

    ! tie-gcm
    use cons_module, only: pi

    implicit none

    ! args
    real, dimension(:,:), intent(inout) :: ocoord_p
    type(reg_grid_dataset_group), intent(inout) :: obs

    ! local
    real(ESMF_KIND_R8), dimension(:), pointer :: lon, lat, alt
    integer :: i, j
    integer :: rc
    integer :: ilon, ilat, ialt

    call ESMF_GridGetCoord( obs%root%grid, coordDim=1, localDE=0, &
        farrayPtr=lon, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    call ESMF_GridGetCoord( obs%root%grid, coordDim=2, localDE=0, &
        farrayPtr=lat, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    call ESMF_GridGetCoord( obs%root%grid, coordDim=3, localDE=0, &
        farrayPtr=alt, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)

    ! observation vector is constructed by reshaping ESMF_field to a flat array
    ! ESMF field has dimension LON x LAT x ALT
    !  vector is flattened using reshape with default order (first dimension varies fastest)
    i = 1
    do j = 1, obs%bundle%size()
      do ialt = lbound(alt,dim=1), ubound(alt,dim=1)
        do ilat = lbound(lat,dim=1), ubound(lat,dim=1)
          do ilon = lbound(lon,dim=1), ubound(lon,dim=1)
            ocoord_p(1,i) = lon(ilon)/180.*pi
            ocoord_p(2,i) = lat(ilat)/180.*pi
            if(size(ocoord_p,dim=1) >= 3) ocoord_p(3,i) = alt(ialt)
            i = i+1
          end do
        end do
      end do
    end do
  end subroutine get_ocoord_lon_lat_lev

  !> Interpolates each observation grid point's geometric position to continuous TIE-GCM grid-cell-index coordinates.
  subroutine get_ocoord_tiegcm_grid_cell_idx_lin(zg,zglevel,ocoord_p,obs,verbose)
    ! computes the position of observation w.r.t to TIE-GCM grid id starting at (0,0,0)
    !
    !   |--1--|--2--|--3--|--4--|  cell index
    !   0     1     2     3     4  coordinates (cell boundaries)
    !   | 0.5 | 1.5 | 2.5 | 3.5 |  cell center coordinate
    !
    ! for periodic dimension
    !   |--1--|--2--|--3--|--4--|--1--|  cell index
    !   0     1     2     3     4     0  coordinates (cell boundaries)
    !   | 0.5 | 1.5 | 2.5 | 3.5 | 0.5 |  cell center coordinate


    ! TODO currently coordinates outside of TIE-GCM mesh are extrpolated.
    ! Instead One could compute ZG for more presure levels!

      ! extern
      use ESMF

      ! tie-gcm
      use cons_module, only: pi
      use fields_module,only: levd0,levd1,lond0,lond1,latd0,latd1

      ! intern
      use array_print_module, only: printMat
      use array_mapping_module, only: flatten
      use cell_id_coordinate_system, only: cell_id_meshgrid
      use quantity_info_module, only: quantity_info, LEVEL_MID, LEVEL_INT, STEP_CURRENT
      use tiegcm_optimized_interpolator, only: state_interpolator

      implicit none

      ! arguments
      real, dimension(:,:,:), intent(inout) :: zg ! size as in fields module
      integer, intent(in) :: zglevel ! LEVEL_MID or LEVEL_INT
      real, dimension(:,:), intent(inout) :: ocoord_p
      type(reg_grid_dataset_group), intent(inout) :: obs
      logical, intent(in), optional :: verbose

      ! local
      type(state_interpolator) :: interpolator
      type(quantity_info) :: info

      real(ESMF_KIND_R8), dimension(:,:,:), pointer, contiguous :: dst_ptr
      type(ESMF_Field) :: dst_field

      integer :: rc
      integer :: dim_id,i,j

      integer :: a,b

      real, dimension(3,size(ocoord_p,dim=2)) :: ocoord_p_geom ! only for control purpose

      real(ESMF_KIND_R8), dimension(levd0:levd1,lond0:lond1,latd0:latd1) :: M

      dst_field  = ESMF_FieldCreate( obs%root%grid, &
                                      typekind = ESMF_TYPEKIND_R8, &
                                          name = "cell_idx",&
                                          rc=rc)

      call ESMF_FieldGet(field=dst_field, localDe=0, farrayPtr=dst_ptr)
      if(ESMF_LogFoundError(rc,msg="get_ocoord_tiegcm_grid_cell_idx:ESMF_FieldGet", &
          rcToReturn=rc)) then
              call shutdown('ESMF error in get_ocoord_tiegcm_grid_cell_idx:ESMF_FieldGet')
      end if

      select case(zglevel)
      case(LEVEL_MID)
        call interpolator%init(dst_field, &
                        zg_mid=zg, &
                        degree=1)
      case(LEVEL_INT)
        call interpolator%init(dst_field, &
                zg_int=zg, &
                degree=1)
      end select

      ! interpolate cell center id for each dimension ----------------------------
      ! --------------------------------------------------------------------------

      do dim_id=1,3

        call cell_id_meshgrid(M,dim_id)

        info%name = "idx"
        info%unit = ""
        info%level=zglevel
        info%step =STEP_CURRENT
        call interpolator%interpolate(info,M,dst_field)

        a = 1
        b = obs%root%n
        do j = 1, obs%bundle%size()
          ocoord_p(dim_id,a:b) = reshape(dst_ptr,(/size(dst_ptr)/))
          a=b+1
          b=j*obs%root%n
        end do
      end do

      call get_ocoord_lon_lat_lev(ocoord_p_geom,obs)
      ocoord_p_geom(1:2,:) = ocoord_p_geom(1:2,:)/pi*180

      if(present(verbose)) then
        if(verbose) then
          write(*,*) 'coordinates of observations w.r.t TIE-GCM grid  (lev lon lat) and geometric position (alt lon lat)'
          do i  = lbound(ocoord_p,dim=2), ubound(ocoord_p,dim=2)
            write(*,'(i6,a,3f6.2,a,3f9.1)') i, "|", ocoord_p(:,i),  "|", ocoord_p_geom(3,i), ocoord_p_geom(1:2,i)
          end do
        end if
      end if

      ! clean --------------------------------------------------------------------
      ! --------------------------------------------------------------------------
      call interpolator%deallocate()
      call ESMF_FieldDestroy(dst_field,noGarbage=.true.)

  end subroutine get_ocoord_tiegcm_grid_cell_idx_lin

end module grid_observation_module
