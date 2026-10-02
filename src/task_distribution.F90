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
! Small helper routines to split a workload count evenly across parallel tasks/threads and compute offsets.

module task_distribution_module
  implicit none

  contains

  !> Splits total_size as evenly as possible across n_threads, returning the
  !! count assigned to each thread.
  subroutine distribute( n_threads, total_size, size_per_thread )

    ! arguments
    integer, intent(in) :: n_threads, total_size
    integer, intent(out), dimension(n_threads) :: size_per_thread

    ! local
    integer :: i
    integer :: npt

    npt = total_size / n_threads

    do i = 1, n_threads
        size_per_thread(i) = npt
        if( i <= total_size - npt*n_threads ) then
            size_per_thread(i) = size_per_thread(i) + 1
        end if
    end do

!         write(*,*) 'size ', total_size, 'distributed on ', n_threads, ' threads : ', size_per_thread

  end subroutine distribute

  !> Computes the starting offset of each thread's chunk in the global array.
  function get_offset( size_per_thread, start_idx ) result( offset )

    ! arguments
    integer, intent(in), dimension(:) :: size_per_thread
    integer, intent(in) :: start_idx

    ! result
    integer, dimension( size(size_per_thread) ) :: offset

    ! local
    integer :: n_threads, i

    n_threads = size(size_per_thread)

    offset(1) = start_idx;
    do i = 2, n_threads
      offset(i) = offset(i-1) + size_per_thread(i-1)
    end do

  end function get_offset

end module task_distribution_module
