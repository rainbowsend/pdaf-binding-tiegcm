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
! Provides field/rank index mapping between the flattened state vector and per-field, per-process array layouts.

! Armin Corbin
! University of Bonn
! APMG
! 26.01.2021
!
! This module (class) provides the mapping between n-dimensional fields,
! that are distributed over several ranks to a 1 dimensional vector
! containing all data and the mapping to a vector containing local data
!
! You can order either by field (F) or rank (R) (see Example below)
!
! order by field        ordered by rank
!
!      | TN  1           | 1 TN
!      | TN  2           | 1 O1
!      | TN  3           | 2 TN
!      | TN  4           | 2 O1
!      | O1  1           | 3 TN
!      | O1  2           | 3 O1
!      | O1  3           | 4 TN
!      | O1  4           | 4 O1
!
! Both methods are useful. F-order is useful for reading hist fields
! and ensemble generation (function 'PDAF_eofcovar'). R-order is useful
! for distributing the state vector with MPI.
!
! ATTENTION R related indices start at 0, but F-related indices start at 1
!
module mapping_module

    IMPLICIT NONE

    type bounds
        integer :: begin    ! start index in global vector
        integer :: back     !   end index in global vector
        integer :: begin_p  ! start index in  local vector
        integer :: back_p   !   end index in  local vector
        integer :: n        ! number of elements between begin and back
    contains
        procedure, pass(this), private :: to_char_array => bounds_to_char_array
    end type bounds

    type mapping

      integer :: ntask    ! number of tasks
      integer :: nfields  ! number of fields
      integer :: n        ! total size of global flattend vector

      integer, dimension(:), allocatable :: size_F    ! global number of elements in field
      integer, dimension(:), allocatable :: size_R    ! number of elements in local vector
      integer, dimension(:,:), allocatable :: F_R_size  ! number of elements in subdomain of field

      !           R1 R2 R3 R4 R5 ...RM    size_F
      !     F1   |__|__|__|__|__|__|__|    |__|
      !     F2   |__|__|__|__|__|__|__|    |__|
      !     ...  |__|__|__|__|__|__|__|    |__|
      !     FN   |__|__|__|__|__|__|__|    |__|
      ! 
      ! size_R   |__|__|__|__|__|__|__|     n

      ! vector starts at 1. first element in offsets is 1
      integer, dimension(:), allocatable :: offsets_F ! start idx of fields in F-ordered vector
      integer, dimension(:), allocatable :: offsets_R ! start idx of fields in R-ordered vector

      ! (1:#fields,0:#task-1)
      type (bounds), dimension(:,:), allocatable :: idx_R ! first order Rank (second field)

      type (bounds), dimension(:,:), allocatable :: idx_F ! first order Field (second rank)

      character(len=16), dimension(:), allocatable :: fd_name
      character(len=64) :: vector_name

      contains
      procedure, pass(this), public :: print => print_mapping
      procedure, pass(this), public :: global_bounds => get_global_bounds
      procedure, pass(this), public :: f2r => rearange_F_2_R
      procedure, pass(this), public :: r2f => rearange_R_2_F
      procedure, pass(this), public :: size => get_size
      procedure, pass(this), public :: deallocate => deconstruct_mapping

    end type mapping

    ! indices used to access state vector cannot be smaller than this value.
    ! usually index starts at one, thus we set it to zero
    integer, parameter :: not_an_index = 0

    contains

    !> Writes the member variables of a bounda instance to a string for printing.
    function bounds_to_char_array(this) result (char_arry)
      implicit none

      class(bounds) :: this

      character(len=60) :: char_arry

      character(len=12), dimension(5) :: out

      write(out(1),'(i12)') this%n

      if(this%begin == not_an_index) then
        write(out(2),'(a12)') '-'
      else
        write(out(2),'(i12)') this%begin
      end if

      if(this%back == not_an_index) then
        write(out(3),'(a12)') '-'
      else
        write(out(3),'(i12)') this%back
      end if

      if(this%begin_p == not_an_index) then
        write(out(4),'(a12)') '-'
      else
        write(out(4),'(i12)') this%begin_p
      end if

      if(this%begin_p == not_an_index) then
        write(out(5),'(a12)') '-'
      else
        write(out(5),'(i12)') this%back_p
      end if

      write(char_arry,'(5a12)') out

    end function

    !> Builds a mapping instance from F_R_size(field, rank) and the associated field names.
    !! The same data is laid out into two differently-ordered flattened vectors of equal total size: idx_R
    !! (rank-major, the layout used to distribute/collect the state vector over MPI)
    !! and idx_F (field-major, the layout used for TIE-GCM history I/O and PDAF_eofcovar
    !! ensemble routines), each with its own begin/back (global) and begin_p/back_p
    !! (local, per rank or per field) bounds.
    function construct_mapping(F_R_size, fd_name, vector_name, offset) result(this)
      !
      ! F_R_size: Matrix with number of elements per field and rank
      !
      !           R1 R2 R3 R4 R5 ...RM
      !     F1   |__|__|__|__|__|__|__|
      !     F2   |__|__|__|__|__|__|__|
      !     ...  |__|__|__|__|__|__|__|
      !     FN   |__|__|__|__|__|__|__|

      implicit none

      ! arguments
      integer, dimension(:,0:), intent(in) :: F_R_size
      character(len=*), dimension(:), intent(in) :: fd_name
      character(len=*), intent(in), optional :: vector_name
      integer, intent(in), optional :: offset

      ! returns
      type(mapping) :: this

      ! local
      integer :: i, gc, lc, tid
      integer :: o0, o1



      if( present( vector_name ) ) then
        write( this%vector_name,'(a)') trim(vector_name)
      else
        this%vector_name=''
      end if

      this%nfields = size(F_R_size, 1)
      this%ntask = size(F_R_size, 2)

      allocate(this%fd_name(this%nfields))
      do i=1, this%nfields
        write(this%fd_name(i),'(a)') trim(fd_name(i))
      end do

      allocate( this%F_R_size(1:this%nfields, 0:this%ntask-1) )
      this%F_R_size(:,:) = F_R_size

      write(*,'(a,i4,a,i4,a)') 'compute mapping for a structure with', this%nfields, ' fields and', this%ntask , ' threads'

      allocate( this%size_F(1:this%nfields) )
      allocate( this%size_R(0:this%ntask-1) )

      allocate( this%idx_R(1:this%nfields, 0:this%ntask-1) )

      ! set global counter
      if(present(offset))then
        if(offset<not_an_index) stop 'invalid offset'
        gc = offset
      else
        gc = 0
      end if

      do tid = 0, this%ntask-1
          lc = 0 ! local counter
          do i = 1, this%nfields
              this%idx_R( i, tid )%n = F_R_size( i, tid )

              if(this%idx_R( i, tid )%n > 0) then

                this%idx_R( i, tid )%begin   = gc + 1
                this%idx_R( i, tid )%begin_p = lc + 1

                gc = gc + this%idx_R( i, tid )%n
                lc = lc + this%idx_R( i, tid )%n

                this%idx_R( i, tid )%back   = gc
                this%idx_R( i, tid )%back_p = lc
              else
                this%idx_R( i, tid )%begin = not_an_index
                this%idx_R( i, tid )%begin_p = not_an_index
                this%idx_R( i, tid )%back = not_an_index
                this%idx_R( i, tid )%back_p = not_an_index
              end if

          end do
          this%size_R( tid ) = lc
      end do

      allocate( this%idx_F(1:this%nfields, 0:this%ntask-1) )

      ! set global counter
      if(present(offset))then
        if(offset<not_an_index) stop 'invalid offset'
        gc = offset
      else
        gc = 0
      end if

      do i = 1, this%nfields
          lc = 0 ! local counter
          do tid = 0, this%ntask-1
              this%idx_F( i, tid )%n = F_R_size( i, tid )

            if( this%idx_F( i, tid )%n > 0) then
              this%idx_F( i, tid )%begin   = gc + 1
              this%idx_F( i, tid )%begin_p = lc + 1

              gc = gc + this%idx_F( i, tid )%n
              lc = lc + this%idx_F( i, tid )%n

              this%idx_F( i, tid )%back   = gc
              this%idx_F( i, tid )%back_p = lc
            else
              this%idx_F( i, tid )%begin = not_an_index
              this%idx_F( i, tid )%begin_p = not_an_index
              this%idx_F( i, tid )%back = not_an_index
              this%idx_F( i, tid )%back_p = not_an_index
            end if
          end do
          this%size_F( i ) = lc
      end do

      call this%global_bounds(begin=o0, back=this%n, fid=this%nfields)

      allocate(this%offsets_F(1:this%nfields))
      do i=1, this%nfields
          call this%global_bounds(begin=this%offsets_F(i), back=o1, fid=i)
      end do

      allocate(this%offsets_R(0: this%ntask-1))
      do i=0, this%ntask-1
          call this%global_bounds(begin=this%offsets_R(i), back=o1, rank=i)
      end do

    end function construct_mapping

    !> Deallocates all arrays owned by this mapping instance.
    subroutine deconstruct_mapping( this )

      implicit none
      class(mapping), intent(inout) :: this

      if (allocated (this%size_F))   deallocate (this%size_F)
      if (allocated (this%size_R))   deallocate (this%size_R)
      if (allocated (this%offsets_F))   deallocate (this%offsets_F)
      if (allocated (this%offsets_R))   deallocate (this%offsets_R)
      if (allocated (this%idx_R))   deallocate (this%idx_R)
      if (allocated (this%idx_F))   deallocate (this%idx_F)
      if (allocated (this%fd_name))   deallocate (this%fd_name)
      if (allocated (this%F_R_size))   deallocate (this%F_R_size)

    end subroutine deconstruct_mapping

    !> Returns the element count for a given rank, a given field, or their intersection, depending on which optional argument is present.
    function get_size(this, rank, fid) result(n)
        implicit none
        ! args
        class(mapping), intent(in) :: this
        integer, intent(in), optional :: rank
        integer, intent(in), optional :: fid

        ! result
        integer :: n

        if(present(rank).and.present(fid))then
          n = this%F_R_size(fid,rank)
        else if (present(fid))then
          n = this%size_F(fid)
        else if (present(rank))then
          n = this%size_R(rank)
        else
          write(*,*) 'ERROR inavlid call of get_size'
          n=0
        end if
    end function


    !> Returns the global start/end indices spanning all fields of a given rank, or all ranks of a given field.
    subroutine get_global_bounds(this, begin, back, rank, fid)

        implicit none

        ! args
        class(mapping), intent(in) :: this
        integer, intent(out) :: begin
        integer, intent(out) :: back
        integer, intent(in), optional :: rank
        integer, intent(in), optional :: fid

        ! local
        character(len=256) :: ermsg
        integer :: i

        ! TODO no error when rank and fid are given

        if(present(rank)) then
          if((rank > this%ntask-1) .or. (rank < 0)) then
            write(ermsg,'(a,i3,a)') 'rank ', rank, ' out of bounds.'
            call shutdown(ermsg)
          end if

          begin = not_an_index
          back  = not_an_index

          do i =1,this%nfields
            if( this%idx_R(i,rank)%begin/=not_an_index) then
              begin = this%idx_R(i,rank)%begin
              exit
            end if
          end do
          do i =this%nfields,1,-1
            if( this%idx_R(i,rank)%back/=not_an_index ) then
              back = this%idx_R(i,rank)%back
              exit
            end if
          end do

        else if (present(fid)) then
          if((fid > this%nfields) .or. (fid < 1)) then
            write(ermsg,'(a,i3,a)') 'index ', fid, ' out of bounds.'
            call shutdown(ermsg)
          end if
          do i =0,this%ntask-1
            if( this%idx_F(fid,i)%begin/=not_an_index ) then
              begin = this%idx_F(fid,i)%begin
              exit
            end if
          end do
          do i =this%ntask-1,0,-1
            if( this%idx_F(fid,i)%back/=not_an_index ) then
              back = this%idx_F(fid,i)%back
              exit
            end if
          end do

        else
          call shutdown('invalid options. either use fid or rank')
        end if

    end subroutine get_global_bounds

    !> Rearranges a global state vector in-place from field-order (F) to rank-order (R) layout.
    subroutine rearange_F_2_R(this, state)
        ! rearanges the the state vector form F to R-order
        ! operates on global state vector

        implicit none

        ! args
        class(mapping), intent(in) :: this
        real, dimension( this%n ), intent(inout) :: state

        ! local
        real, dimension( this%n  ) :: state_copy
        integer :: i, tid

        state_copy = state

        do tid = 0, this%ntask-1
            do i = 1, this%nfields
                if(this%idx_R(i,tid)%n > 0) then
                       state( this%idx_R(i,tid)%begin : this%idx_R(i,tid)%back ) = &
                  state_copy( this%idx_F(i,tid)%begin : this%idx_F(i,tid)%back )
                end if
            end do
        end do

    end subroutine rearange_F_2_R

    !> Rearranges a global state vector in-place from rank-order (R) to field-order (F) layout.
    subroutine rearange_R_2_F(this, state)
        ! rearanges the the state vector form R to F-order
        ! operates on global state vector

        implicit none

        ! args
        class(mapping), intent(in) :: this
        real, dimension( this%n ), intent(inout) :: state

        ! local
        real, dimension( this%n ) :: state_copy
        integer :: i, tid

        state_copy = state

        do tid = 0, this%ntask-1
            do i = 1, this%nfields
                if(this%idx_F(i,tid)%n > 0) then
                       state( this%idx_F(i,tid)%begin : this%idx_F(i,tid)%back ) = &
                  state_copy( this%idx_R(i,tid)%begin : this%idx_R(i,tid)%back )
                end if
            end do
        end do

    end subroutine rearange_R_2_F

    !> Prints the full rank-order and field-order vector layouts (per-field/per-rank bounds, offsets, sizes) to stdout for debugging.
    subroutine print_mapping( this )

      implicit none

      class(mapping), intent(in) :: this

      integer :: tid, i
      integer :: o0, o1

      write(*, '(2a)') trim(this%vector_name), ' VECTOR LAYOUT -- RANK (SUBDOMAIN) FIRST'
      write(*, '(2a)') repeat(" ", 33), '|          GLOBAL       |          LOCAL'
      write(*, '(a)') '           field   task         n       start         end' &
                    //'       start         end'
      do tid = 0, this%ntask-1
          call this%global_bounds(begin=o0, back=o1,rank=tid)
          write(*,'(a,3i12)')  repeat(" ", 21), this%size_R(tid), o0, o1
          do i = 1, this%nfields
          write(*,'(a16, i5, a60)') &
                        trim(this%fd_name(i)), tid, &
                        this%idx_R( i, tid )%to_char_array()
          end do
          write(*,'(a)') repeat("-", 81)
      end do
      write(*,'(a)') repeat("-", 81)
      write(*, '(a)') '  rank start/offset   count/size'
      do tid = 0, this%ntask-1
        write(*,'(i6,2i13)') tid, this%offsets_R(tid), this%size_R(tid)
      end do

      write(*, '(2a)') trim(this%vector_name), ' VECTOR LAYOUT -- FIELD FIRST'
      write(*, '(2a)') repeat(" ", 33), '|          GLOBAL       |          LOCAL'
      write(*, '(a)') '           field  task          n       start         end' &
                    //'       start         end'
      do i = 1, this%nfields
          call this%global_bounds(begin=o0, back=o1,fid=i)
          write(*,'(a,3i12)') repeat(" ", 21), this%size_F(i), o0, o1
          do tid = 0, this%ntask-1
          write(*,'(a16, i5, a60)') &
                        trim(this%fd_name(i)),  tid, &
                         this%idx_F( i, tid )%to_char_array()
          end do
          write(*,'(a)') repeat("-", 81)
      end do
      write(*,'(a)') repeat("-", 81)
      write(*, '(a)') '           field start/offset   count/size'
      do i = 1, this%nfields
        write(*,'(a16,2i13)') trim(this%fd_name(i)), this%offsets_F(i), this%size_F(i)
      end do

    end subroutine print_mapping

end module mapping_module
