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
! Generic named-field container type ("bundle") with typed scalar/1D/2D/3D accessors for grouped state/observation data.

module field_bundle_module

use uset_module, only: char_uset

implicit none


  type :: bundle
    type(char_uset), private :: names
    real, dimension(:,:), allocatable :: data
    integer, allocatable, dimension(:) :: shape
    contains
    procedure, pass(this) :: set_shape => bundle_set_shape
    procedure, pass(this) :: allocate => bundle_allocate
    procedure, pass(this) :: deallocate => bundle_deallocate
    procedure, pass(this) :: size => bundle_size
    procedure, pass(this) :: get_name => bundle_get_name
    procedure, pass(this) :: get_id => bundle_get_id
    procedure, pass(this) :: bundle_get_scalar
    procedure, pass(this) :: bundle_get_1d
    procedure, pass(this) :: bundle_get_2d
    procedure, pass(this) :: bundle_get_3d
    generic, public :: get => bundle_get_scalar, &
                              bundle_get_1d, &
                              bundle_get_2d, &
                              bundle_get_3d
    procedure, pass(this) :: get_flat => bundle_get_flat
  end type

  contains

  !> Reassigns a pointer to an existing 2D array without copying (used to obtain a
  !! pointer to the bundle's underlying data).
  !!
  !! TODO: is this an evil hack?
  subroutine get_ptr_2d(in,out)

    implicit none
    real, dimension(:,:), contiguous, target, intent(in) :: in
    real, dimension(:,:), contiguous, pointer, intent(out) :: out

    out => in
  end subroutine

  !> Sets (or, if omitted, defaults to scalar) the per-field data shape used when
  !! allocating this bundle.
  subroutine bundle_set_shape(this,data_shape)
    class(bundle) :: this
    integer, dimension(:), intent(in), optional :: data_shape

    if(present(data_shape)) then
      allocate(this%shape, source=data_shape)
    else
      ! assume scalar
      allocate(this%shape(1))
      this%shape(1) = 1
    end if
  end subroutine

  !> Allocates the bundle's flat data array for the given field names (and shape),
  !! initialized to TIE-GCM's fill/missing value.
  subroutine bundle_allocate(this,field_names,data_shape)

    use params_module,only: spval

    implicit none


    class(bundle) :: this
    type(char_uset), intent(in) :: field_names
    integer, dimension(:), intent(in), optional :: data_shape

    if(present(data_shape)) then
      call this%set_shape(data_shape)
    end if

    if(allocated(this%shape).eqv..false.)then
      call shutdown('you need to set the shape of the bundle before allocating it')
    end if

    if(allocated(this%data).eqv..true.)then
      call shutdown('bunde has already been allocated')
    end if

    this%names = field_names

    write(*,*) "allocating field bundle of size", this%names%size(), &
               "and local shape ", this%shape

    allocate(this%data(product(this%shape,dim=1),this%names%size()))
    this%data = spval

  end subroutine

  !> Deallocates the bundle's field names, data array, and shape.
  subroutine bundle_deallocate(this)

    implicit none
    class(bundle) :: this

    call this%names%deallocate
    if(allocated(this%data)) deallocate(this%data)
    if(allocated(this%shape)) deallocate(this%shape)

  end subroutine


  !> Returns the number of fields stored in the bundle.
  function bundle_size(this) result(res)
    implicit none

    ! arguments
    class(bundle) :: this

    ! result
    integer :: res

    res = this%names%size()

  end function

  !> Returns the index of a named field within the bundle, aborting if not found.
  function bundle_get_id(this, field_name) result(id)

    implicit none

    ! arguments
    class(bundle) :: this
    character(len=*), intent(in) :: field_name

    ! result
    integer :: id

    id = this%names%findloc(field_name)
    if(id<=0)then
      call shutdown(trim(field_name)//' is not in field_bundle')
    end if

  end function

  !> Returns the name of the field at a given index in the bundle.
  function bundle_get_name(this,idx) result(res)
    implicit none

    ! arguments
    class(bundle) :: this
    integer, intent(in) :: idx
    ! result
    character(len=16) :: res

    res = this%names%at(idx)

  end function

 !> Returns a flat pointer to a named field's data within the bundle.
 subroutine bundle_get_flat(this, field_name, flat)
    implicit none

    ! arguments
    class(bundle) :: this
    character(len=*), intent(in) :: field_name
    real, dimension(:), contiguous, pointer, intent(out) :: flat

    ! local
    integer :: id
    real, dimension(:,:), contiguous, pointer :: bundle_data


    if(allocated(this%data).eqv..false.)then
      call shutdown('tried to access unallocated bundle')
    end if

    id = this%get_id(field_name)
    call get_ptr_2d(this%data,bundle_data)

    flat => bundle_data(:,id)

  end subroutine

  !> Returns a scalar pointer to a named field; requires the bundle's shape to be scalar.
  subroutine bundle_get_scalar(this, field_name, field_ptr)
    ! arguments
    class(bundle) :: this
    character(len=*), intent(in) :: field_name
    real, pointer, intent(out) :: field_ptr

    ! local
    real, dimension(:), contiguous, pointer :: ptr

    if(product(this%shape,dim=1)/=1 )then
      call shutdown("Requested a scalar, but field bundle is not scalar")
    end if

    call this%get_flat(field_name,ptr)
    field_ptr => ptr(1)

  end subroutine

  !> Returns a 1D-reshaped pointer to a named field; requires the bundle's shape to
  !! be rank 1.
  subroutine bundle_get_1d(this, field_name, field_ptr)
    implicit none

    ! arguments
    class(bundle) :: this
    character(len=*), intent(in) :: field_name
    real, dimension(:), contiguous, pointer, intent(out) :: field_ptr


    ! local
    real, dimension(:), contiguous, pointer :: ptr

    if(size(this%shape,dim=1)/=1) call shutdown('rank mismatch expected rank 1')

     call this%get_flat(field_name,ptr)
     field_ptr(1:this%shape(1)) => ptr

  end subroutine

  !> Returns a 2D-reshaped pointer to a named field; requires the bundle's shape to
  !! be rank 2.
  subroutine bundle_get_2d(this, field_name, field_ptr)
    implicit none

    ! arguments
    class(bundle) :: this
    character(len=*), intent(in) :: field_name
    real, dimension(:,:), contiguous, pointer, intent(out) :: field_ptr


    ! local
    real, dimension(:), contiguous, pointer :: ptr

    if(size(this%shape,dim=1)/=2) call shutdown('rank mismatch expected rank 2')

     call this%get_flat(field_name,ptr)
     field_ptr(1:this%shape(1),1:this%shape(2)) => ptr

  end subroutine

  !> Returns a 3D-reshaped pointer to a named field; requires the bundle's shape to
  !! be rank 3.
  subroutine bundle_get_3d(this, field_name, field_ptr)
    implicit none

    ! arguments
    class(bundle) :: this
    character(len=*), intent(in) :: field_name
    real, dimension(:,:,:), contiguous, pointer, intent(out) :: field_ptr


    ! local
    real, dimension(:), contiguous, pointer :: ptr

    if(size(this%shape,dim=1)/=3) call shutdown('rank mismatch expected rank 3')

     call this%get_flat(field_name,ptr)
     field_ptr(1:this%shape(1),1:this%shape(2),1:this%shape(3)) => ptr

  end subroutine

end module
