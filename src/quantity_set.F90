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
! Bookkeeping for named sets of quantities (state, observed, mandatory, prognostic, diagnostic, output) used to decide what is assimilated/written.

module quantity_sets_module

use uset_module, only: char_uset

implicit none

type :: quantity_sets
  type(char_uset) :: mandatory
  type(char_uset) :: state
  type(char_uset) :: observed
  ! the following sets are derived from the above
  type(char_uset) :: total
  type(char_uset) :: prognostic
  type(char_uset) :: diagnostic
  type(char_uset) :: output
  contains
  procedure, pass(this) :: deallocate => quantity_sets_deallocate
  procedure, pass(this) :: add => quantity_sets_add
  procedure, pass(this) :: update_output => quantity_sets_update_output
  procedure, pass(this) :: print => quantity_sets_print
end type

contains

!> Deallocates all member variables of this instance of a quantity_sets.
subroutine quantity_sets_deallocate(this)

  implicit none

  ! arguments
  class(quantity_sets) :: this

  call this%mandatory%deallocate
  call this%state%deallocate
  call this%observed%deallocate
  call this%total%deallocate
  call this%prognostic%deallocate
  call this%diagnostic%deallocate
  call this%output%deallocate
end subroutine

!> Adds fields to the state/observed/mandatory set, then updates the
!! derived total, prognostic, diagnostic, and output sets.
subroutine quantity_sets_add(this, fields, set)

    ! tie-gcm
    use fields_module,only: f4d

    implicit none

    ! arguments
    class(quantity_sets) :: this
    character(len=*), dimension(:), intent(in) :: fields
    character(len=*), intent(in) :: set

    ! local
    integer :: i
    logical, dimension(size(fields)) :: mask

    if(size(fields)<1)then
      write(*,*) 'Error in update_sets. Input is empty'
      return
    end if

    ! mask empty strings
    do i = 1, size(fields,dim=1)
      mask(i) = fields(i)/= ""
    end do

    select case (set)
      case ("state")
        call this%state%add(PACK(fields,mask))
      case ("observed")
        call this%observed%add(PACK(fields,mask))
      case ("mandatory")
        call this%mandatory%add(PACK(fields,mask))
      case default
        call shutdown("quantity_sets_add: "//set//" is an invalid option")
    end select

    call this%total%deallocate
    call this%total%add((/this%state,this%observed,this%mandatory/))

     do i=1, this%total%size()
      if(findloc( f4d%short_name, trim(this%total%at(i)),dim=1 ) > 0) then
        call this%prognostic%add(this%total%at(i))
      end if
    end do

    this%diagnostic=this%total%difference(this%prognostic)

    call this%update_output()

end subroutine

!> Rebuilds the output set from the state, observed, and mandatory sets
!! according to the cfg_output save flags.
subroutine quantity_sets_update_output(this)

  ! intern
  use character_routines_module, only: string_ends_with
  use configuration, only: cfg_output

  implicit none

  ! arguments
  class(quantity_sets) :: this

  ! local
  integer :: i

  call this%output%deallocate

  if (cfg_output%save_state) then
    if(cfg_output%save_state_nm.eqv..false.) then
      do i = 1, this%state%size()
        if(string_ends_with(trim(this%state%at(i)),"_NM").eqv..false.) then
          call this%output%add(this%state%at(i))
        end if
      end do
    else
      call this%output%add(this%state)
    end if
  end if

  if(cfg_output%save_obs) then
    call this%output%add(this%observed)
  end if

  call this%output%add(this%mandatory)

end subroutine

!> Prints all quantity sets to screen.
subroutine quantity_sets_print(this)

  implicit none

  ! arguments
  class(quantity_sets) :: this

  write(*,'(a,/,3x)',advance="no") 'all quantities:'
  call  this%total%print()

  write(*,'(a,/,3x)',advance="no") 'observed quantities:'
  call  this%observed%print()

  write(*,'(a,/,3x)',advance="no") 'mandatory quantities (specified in configuration file):'
  call  this%mandatory%print()

  write(*,'(a,/,3x)',advance="no") 'diagnostic quantities:'
  call  this%diagnostic%print()

  write(*,'(a,/,3x)',advance="no") 'prognostic quantities:'
  call  this%prognostic%print()

  write(*,'(a,/,3x)',advance="no") 'quantities written to result file:'
  call  this%output%print()

  write(*,'(/)',advance="no")
end subroutine

end module
