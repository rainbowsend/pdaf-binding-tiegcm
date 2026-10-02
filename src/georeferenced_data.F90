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
! Manages point and regular-grid observation datasets, interpolating TIE-GCM state onto them for output/comparison.

module georeferenced_data_module

! intern
use trajectory_data_module, only: trajectory_data
use field_bundle_module, only: bundle
use grid_observation_module, only: reg_grid_dataset_root
use result_file_writer_module, only: nc_trajectory, nc_reg_grid

implicit none

type, abstract :: georeferenced_data_type
  character(len=16) :: name
  type(bundle) :: data
  contains
  procedure(interpolate_interface), deferred :: interpolate
  procedure(deallocate_interface), deferred :: deallocate
end type

abstract interface
  !> Deferred interface: interpolates the given quantities onto this dataset's
  !! points/cells.
  subroutine interpolate_interface(this, quantities, zg_mid, zg_mid_nm, zg_int, zg_int_nm)
    use quantity_info_module, only: quantity_type
    import georeferenced_data_type
    class(georeferenced_data_type) :: this
    type(quantity_type), dimension(:), intent(in) :: quantities
    real, dimension(:,:,:), intent(inout) ::  zg_mid, zg_mid_nm, zg_int, zg_int_nm
  end subroutine
  !> Deferred interface: deallocates this dataset instance.
  subroutine deallocate_interface(this)
    use quantity_info_module, only: quantity_type
    import georeferenced_data_type
    class(georeferenced_data_type) :: this
  end subroutine
end interface

type :: georeferenced_data
  class(georeferenced_data_type), pointer :: ptr => null()
end type

type, extends(georeferenced_data_type) :: point
  type(trajectory_data), pointer :: trajectory => null()
  class(nc_trajectory), pointer :: writer => null()
  contains
  procedure, pass(this) :: init => point_init
  procedure, pass(this) :: interpolate => point_interpolate
  procedure, pass(this) :: deallocate => point_deallocate
end type

type, extends(georeferenced_data_type) :: regular_grid
  type(reg_grid_dataset_root), pointer :: grid => null()
  class(nc_reg_grid), pointer :: writer => null()
  contains
  procedure, pass(this) :: init => regular_grid_init
  procedure, pass(this) :: interpolate => regular_grid_interpolate
  procedure, pass(this) :: deallocate => regular_grid_deallocate
end type

integer, parameter :: max_datasets = 10
integer, protected :: n_datasets = 0
type(georeferenced_data), dimension(max_datasets) :: georeferenced_datasets

contains

!> Registers a georeferenced_data_type instance in the global list of datasets.
subroutine link_georeferenced_data(data)

    implicit none

    ! arguments
    class(georeferenced_data_type), target :: data

    n_datasets = n_datasets + 1
    georeferenced_datasets(n_datasets)%ptr => data

end subroutine

!> Initializes a point-type dataset, linking it to a trajectory and NetCDF writer.
subroutine point_init(this, trajectory, writer, name)

  use result_file_writer_module, only: nc_trajectory

  implicit none

  ! arguments
  class(point) :: this
  type(trajectory_data), intent(in), target :: trajectory
  type(nc_trajectory), intent(in), target :: writer
  character(len=*), intent(in), optional :: name

  if(present(name))then
    write(this%name,*) trim(name)
  end if

  call this%data%set_shape()
  call link_georeferenced_data(this)
  this%trajectory => trajectory
  this%writer => writer

end subroutine

!> Deallocates a point dataset's data bundle and detaches its trajectory/writer.
subroutine point_deallocate(this)

  implicit none

  ! arguments
  class(point) :: this

  call this%data%deallocate()
  ! TODO unlink georeferenced_datasets
  this%trajectory => null()
  this%writer => null()
end subroutine

!> Interpolates the requested quantities onto the trajectory's current along-track
!! position, if the trajectory has valid data at the current model time.
subroutine point_interpolate(this, quantities, zg_mid, zg_mid_nm, zg_int, zg_int_nm)

  ! intern
  use configuration, only: cfg_filter
  use quantity_info_module, only: quantity_type, get_required_zg, supports_point_interpolation
  use tiegcm_optimized_interpolator, only: sparse_state_interpolator
  use time_module, only: get_current_modeltime

  implicit none

  ! arguments
  class(point) :: this
  type(quantity_type), dimension(:), intent(in) :: quantities
  real, dimension(:,:,:), intent(inout) ::  zg_mid, zg_mid_nm, zg_int, zg_int_nm

  ! local
  type(sparse_state_interpolator) :: interp
  logical, dimension(4) :: req_ZG
  real, pointer :: dst
  real :: val
  real :: modeltime
  logical :: is_available
  integer :: i
  real, dimension(3) :: position

  call get_current_modeltime(modeltime)
  is_available = this%trajectory%valid_at_epoch(modeltime)

  if(is_available)then

    call get_required_zg(quantities,&
      req_mid=req_ZG(1),&
      req_mid_nm=req_ZG(2),&
      req_int=req_ZG(3),&
      req_int_nm=req_ZG(4))

    call this%trajectory%get_at_epoch(modeltime,val,position)

    call interp%init(position,&
                    zg_mid=zg_mid,&
                    zg_mid_nm=zg_mid_nm,&
                    zg_int=zg_int,&
                    zg_int_nm=zg_int_nm,&
                    init_flags=req_ZG,&
                    degree=cfg_filter%spline_degree )

    do i=1,size(quantities)

      if(supports_point_interpolation(quantities(i)%info%name).eqv..false.) cycle

      dst => null()
      call this%data%get(quantities(i)%info%name, dst)

!       if(associated(dst).eqv..false.) then
!         call shutdown("point_interpolate: invalid pointer for "// quantities(i)%info%name)
!       end if

      call interp%interpolate(quantities(i)%info,&
                              quantities(i)%info%data,&
                              dst)
    end do

    call interp%deallocate()
  end if
end subroutine

!> Initializes a regular-grid dataset, linking it to a grid and NetCDF writer.
subroutine regular_grid_init(this, grid, writer, name)

  use result_file_writer_module, only: nc_reg_grid

  implicit none

  ! arguments
  class(regular_grid) :: this
  type(reg_grid_dataset_root), intent(in), target :: grid
  type(nc_reg_grid), intent(in), target :: writer
  character(len=*), intent(in), optional :: name

  ! local
  integer, dimension(3) :: grid_shape

  if(present(name))then
    write(this%name,*) trim(name)
  end if

  ! TODO move grid_shape init to function in grid_obs
  grid_shape(1)=grid%nlon
  grid_shape(2)=grid%nlat
  grid_shape(3)=grid%nalt

  call this%data%set_shape(grid_shape)
  call link_georeferenced_data(this)
  this%grid => grid
  this%writer => writer


end subroutine

!> Deallocates a regular-grid dataset's data bundle and detaches its grid/writer.
subroutine regular_grid_deallocate(this)

  implicit none

  ! arguments
  class(regular_grid) :: this

  call this%data%deallocate()
  ! TODO unlink georeferenced_datasets
  this%grid => null()
  this%writer => null()
end subroutine

!> Interpolates the requested quantities onto the regular output grid via ESMF fields.
subroutine regular_grid_interpolate(this, quantities, zg_mid, zg_mid_nm, zg_int, zg_int_nm)

  ! extern
  use esmf

  ! intern
  use configuration, only: cfg_filter
  use quantity_info_module, only: quantity_type, get_required_zg, quantity_type
  use tiegcm_optimized_interpolator, only: state_interpolator

  implicit none

  ! arguments
  class(regular_grid) :: this
  type(quantity_type), dimension(:), intent(in) :: quantities
  real, dimension(:,:,:), intent(inout) ::  zg_mid, zg_mid_nm, zg_int, zg_int_nm

  ! local
  type(ESMF_Field) :: dst_field
  real(ESMF_KIND_R8), dimension(:,:,:), contiguous, pointer :: dst_ptr
  type(state_interpolator) interp

  integer :: i
  integer :: rc

  logical, dimension(4) :: req_ZG

  call get_required_zg(quantities,&
    req_mid=req_ZG(1),&
    req_mid_nm=req_ZG(2),&
    req_int=req_ZG(3),&
    req_int_nm=req_ZG(4))

  call interp%init( &
            grid=this%grid%grid, &
            zg_mid=zg_mid,&
            zg_mid_nm=zg_mid_nm,&
            zg_int=zg_int,&
            zg_int_nm=zg_int_nm,&
            init_flags=req_ZG,&
            degree=cfg_filter%spline_degree)

  do i=1,size(quantities)

    dst_ptr => null()
    call this%data%get(quantities(i)%info%name, dst_ptr)

!     if(associated(dst_ptr).eqv..false.) then
!       call shutdown("regular_grid_interpolate: invalid pointer for "// quantities(i)%info%name)
!     end if

    dst_field = ESMF_FieldCreate( this%grid%grid, &
                                  dst_ptr, &
                                  name = "result computatiton dst", &
                                  rc=rc)
    if(ESMF_LogFoundError(rc,msg="regular_grid_interpolate:ESMF_FieldCreate", &
        rcToReturn=rc)) then
        call shutdown('ESMF error in regular_grid_interpolate:ESMF_FieldCreate')
    end if
    call interp%interpolate(quantities(i)%info,&
                            quantities(i)%info%data,&
                            dst_field)
    call ESMF_FieldDestroy(dst_field,noGarbage=.true.)

  end do

  call interp%deallocate()

end subroutine

!> Allocates the data bundle for the given field names in every registered
!! georeferenced dataset.
subroutine georeferenced_data_allocate(names)

  use uset_module, only: char_uset

  implicit none

  ! arguments
  type(char_uset), intent(in) :: names

  ! local
  integer :: i

  do i=1,n_datasets
    ! TODO deallocation
    call georeferenced_datasets(i)%ptr%data%allocate(names)
  end do

end subroutine

!> Copies a regridded field bundle into the flat observation state vector, per the
!! field/rank mapping.
subroutine map_bundle_to_vec(grid_obs, regridded, m_state_p)

    ! extern
    use ESMF
    ! tie-gcm
    use mpi_module, only: mytid

    ! intern
    use field_bundle_module, only: bundle
    use grid_observation_module, only: reg_grid_dataset_group

    implicit none

    ! arguments

    type (reg_grid_dataset_group), intent(in) :: grid_obs
    type(bundle), intent(in) :: regridded
    real, dimension(:), intent(out)  :: m_state_p

    ! local
    integer :: i
    real, contiguous, pointer :: field_data(:)
    write(*,*)  'assining observation field to observation vector'
    do i=1, size(grid_obs%map%fd_name)


        write(*,*) grid_obs%map%fd_name(i), ' ', &
                   grid_obs%map%idx_R(i,mytid)%begin_p, &
                   grid_obs%map%idx_R(i,mytid)%back_p, &
                   grid_obs%map%F_R_size(i,mytid)

        field_data => null()
        call regridded%get_flat(grid_obs%map%fd_name(i), field_data)

        m_state_p( grid_obs%map%idx_R(i,mytid)%begin_p : &
                   grid_obs%map%idx_R(i,mytid)%back_p ) &
                  = field_data

    end do
end subroutine

!> Copies a regridded ESMF field bundle into the flat observation state vector, per
!! the field/rank mapping.
!!
!! TODO: superseded by map_bundle_to_vec; retained for legacy callers.
subroutine assign_to_observation_vector(grid_obs, regridded, m_state_p)

    ! extern
    use ESMF
    ! tie-gcm
    use mpi_module, only: mytid

    ! intern
    use grid_observation_module, only: reg_grid_dataset_group, field_bundle


    implicit none

    ! arguments

    type( field_bundle ) :: regridded

    type (reg_grid_dataset_group), intent(inout) :: grid_obs
    real, dimension(:), intent(inout)  :: m_state_p

    ! local
    integer :: i, rc
    integer :: state_field_id
    real(ESMF_KIND_R8), contiguous, pointer :: field_data(:,:,:)
    write(*,*)  'assining observation field to observation vector'
    do i=1, size(grid_obs%map%fd_name)

        state_field_id = findloc( regridded%names, grid_obs%map%fd_name(i), dim=1 )

        if ( state_field_id < 1) then
            call shutdown('msis field '// grid_obs%map%fd_name(i) // ' not found in regridded state field')
        end if

        field_data => null()
        call ESMF_FieldGet(regridded%bundle(state_field_id), localDe=0, farrayPtr=field_data, &
            rc=rc)

        write(*,*) grid_obs%map%fd_name(i), ' ', &
                   grid_obs%map%idx_R(i,mytid)%begin_p, &
                   grid_obs%map%idx_R(i,mytid)%back_p, &
                   grid_obs%map%F_R_size(i,mytid)

        m_state_p( grid_obs%map%idx_R(i,mytid)%begin_p : &
                   grid_obs%map%idx_R(i,mytid)%back_p ) &
                  = RESHAPE( field_data, (/grid_obs%map%F_R_size(i,mytid)/) )

    end do
end subroutine

end module
