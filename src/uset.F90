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
! Generic dynamically-growing unique-string-set container (add/union/difference/print) used throughout the coupling layer.

!
! Impementation of a unique set. When adding items it is ensured dublicates
! are not added to the data structure.
! Memory is automatically allocated in an efficent way.
module uset_module

  implicit none

  integer, parameter :: str_len_char_uset = 16

  type :: char_uset
      character(len=str_len_char_uset), allocatable, dimension(:), private :: items
      integer, private :: n_items = 0
    contains
      procedure, pass(this), private :: char_uset_add_item
      procedure, pass(this), private :: char_uset_add_items
      procedure, pass(this), private :: char_uset_add_other
      procedure, pass(this), private :: char_uset_add_others
      generic, public :: add => char_uset_add_item,&
                                char_uset_add_items,&
                                char_uset_add_other,&
                                char_uset_add_others
      procedure, pass(this), public :: at => char_uset_at
      procedure, pass(this), public :: to_array => char_uset_to_array
      procedure, pass(this), public :: findloc => char_uset_findloc
      procedure, pass(this), public :: has => char_uset_has_item
      procedure, pass(this), public :: size => char_uset_size
      procedure, pass(this), public :: print => char_uset_print
      procedure, pass(this), private :: reallocate => char_uset_reallocate
      procedure, pass(this), public :: deallocate => char_uset_deconstruct
      procedure, pass(this), public :: difference => char_uset_difference
      procedure, pass(this), public :: union => char_uset_union
  end type char_uset

  contains

  !> Returns the number of items currently stored in the set.
  function char_uset_size(this) result(s)
    implicit none
    ! arguments
    class(char_uset), intent(in) :: this
    ! returns
    integer:: s

    s=this%n_items
  end function char_uset_size

  !> Returns the index of item in the set, or 0 if it is not present.
  function char_uset_findloc(this,item) result(loc)

    implicit none

    ! input arguments
    class(char_uset), intent(in) :: this
    character(len=*), intent(in) :: item

    ! local
    integer :: loc

    loc=0
    if (allocated(this%items).and.(size(this%items)>0)) then
        loc = findloc(this%items(1:this%n_items),item,dim=1)
    end if

  end function char_uset_findloc

  !> Returns whether item is already present in the set.
  function char_uset_has_item(this,item) result(res)

    implicit none

    ! input arguments
    class(char_uset), intent(in) :: this
    character(len=*), intent(in) :: item

    ! result
    logical :: res

    ! local
    integer :: loc

    loc = this%findloc(item)

    res = loc > 0

  end function char_uset_has_item

  !> Returns the items in this that are not present in other.
  function char_uset_difference(this,other) result(difference)
    implicit none
    ! arguments
    class(char_uset), intent(inout) :: this
    type(char_uset), intent(in) :: other

    ! result
    type(char_uset) :: difference

    ! local
    integer :: i

    do i = 1, this%n_items
      if( other%has(this%at(i)).eqv..false. )then
        call difference%add(this%at(i))
      end if
    end do

  end function

  !> Returns the union of this and other, without duplicates.
  function char_uset_union(this,other) result(union)
    implicit none
    ! arguments
    class(char_uset), intent(inout) :: this
    type(char_uset), intent(in) :: other

    ! result
    type(char_uset) :: union

    call union%add(other=this)
    call union%add(other=other)

  end function

  !> Adds item to the set, growing storage if needed, unless it is already present.
  subroutine char_uset_add_item(this,item)

    implicit none

    ! arguments
    class(char_uset), intent(inout) :: this
    character(len=*), intent(in) :: item

    if(allocated(this%items).eqv..false.) then
      this%n_items = 0
    end if

    if(char_uset_has_item(this,item).eqv..false.)then
      this%n_items = this%n_items+1
      call char_uset_reallocate(this,this%n_items)
      write(this%items(this%n_items),'(16a)') item
    end if

  end subroutine

  !> Adds each element of items to the set, skipping duplicates.
  subroutine char_uset_add_items(this,items)

    implicit none

    ! arguments
    class(char_uset), intent(inout) :: this
    character(len=*), dimension(:), intent(in) :: items

    ! local
    integer :: i

    if(size(items)> 0) then
      do i=1,size(items,dim=1)
        call char_uset_add_item(this,items(i))
      end do
    end if

  end subroutine

  !> Adds all items of another set to this set, skipping duplicates.
  subroutine char_uset_add_other(this,other)

    implicit none

    ! arguments
    class(char_uset), intent(inout) :: this
    type(char_uset), intent(in) :: other


    if(other%n_items>0)then
      call char_uset_add_items(this,other%items(1:other%n_items))
    end if

  end subroutine

  !> Adds all items of an array of other sets to this set, skipping duplicates.
  subroutine char_uset_add_others(this,others)

    implicit none

    ! arguments
    class(char_uset), intent(inout) :: this
    type(char_uset), dimension(:), intent(in) :: others

    ! local
    integer :: i

    if(size(others)> 0) then
      do i=1,size(others,dim=1)
        call char_uset_add_other(this,others(i))
      end do
    end if

  end subroutine

  !> Returns the item at position idx, printing a warning if idx is out of bounds.
  function char_uset_at(this,idx) result(item)

    implicit none

    ! arguments
    class(char_uset), intent(in) :: this
    integer, intent(in) :: idx

    ! result
    character(len=str_len_char_uset) :: item

    if(idx < 1 .or. idx > this%n_items ) then
      write(*,*) "uset index is out of bounds"
    end if
    item = this%items(idx)

  end function char_uset_at

  !> Copies the set's items into an allocatable array.
  subroutine char_uset_to_array(this,items)

    implicit none

    ! arguments
    class(char_uset), intent(in) :: this
    character(len=str_len_char_uset), allocatable, dimension(:), intent(inout) :: items

    if (allocated(items)) deallocate(items)
    allocate( items, source=this%items(1:this%n_items) )

  end subroutine char_uset_to_array

  !> Deallocates the set's storage and resets it to empty.
  subroutine char_uset_deconstruct(this)
    ! arguments
    class(char_uset), intent(inout) :: this

    if(allocated(this%items)) deallocate(this%items)
    this%n_items=0

  end subroutine char_uset_deconstruct

  !> Writes all items in the set to standard output on one line.
  subroutine char_uset_print(this)
    ! arguments
    class(char_uset), intent(in) :: this

    integer :: i
    ! write(*,*) 'uset(',this%n_items,'): '
    do i=1,this%n_items
      write(*,'(3a,1x)',advance="no") "'",trim(this%items(i)),"'"
    end do
    write(*,'(a)') ''

  end subroutine char_uset_print

  !> Grows the set's internal storage to at least newsize, preserving existing items.
  subroutine char_uset_reallocate(this, newsize)
    ! arguments
    class(char_uset), intent(inout) :: this
    integer, intent(in) :: newsize

    !
    character(len=str_len_char_uset), allocatable, dimension(:) :: tmp

    if(allocated(this%items) .eqv. .false.) then
      allocate(this%items(2*newsize))
    else if(size(this%items) < newsize) then

      allocate( tmp,source=this%items )

      if(allocated(this%items)) deallocate(this%items)
      allocate(this%items(2*newsize))
      this%items(lbound(tmp,dim=1):ubound(tmp,dim=1)) = tmp
      deallocate(tmp)
    end if
  end subroutine char_uset_reallocate

end module uset_module
