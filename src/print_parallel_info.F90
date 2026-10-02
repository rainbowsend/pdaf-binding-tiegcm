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
! Diagnostic printouts of the MPI rank/task layout (world, filter, model, coupling ranks) for debugging parallel setup.

!
MODULE print_parallel_info

  ! extern
  use mpi_f08

  ! tie-gcm
  use mpi_module, only: TIEGCM_WORLD, mytid

  ! intern
  use mod_parallel_pdaf

  implicit none

  contains

  !> Prints this process's world/model/filter/couple ranks and filter-PE status, prefixed by an optional message.
  subroutine print_ranks(msg_in)

    implicit none

    CHARACTER(len=*), intent(in), optional :: msg_in
    CHARACTER(len=1024) :: msg

    if( present( msg_in ) ) then
        msg = msg_in
    else
        msg = ''
    end if

    write(*,'(a/,a,i4/,a,i4/,a,i4/,a,l4/,a,i4)')  msg, &
                            '   world rank: ', rank_world, &
                            ' model run id: ', task_id, &
                            '  rank filter: ', rank_filter, &
                            '     filterPE: ', filterpe, &
                            '  rank couple: ', rank_couple
  end subroutine print_ranks

  !> Gathers each process's world/filter/model/couple ranks and filter-PE flag to rank 0 and prints the full PE configuration table.
  subroutine print_mpi_layout

    implicit none

    integer :: MPIerr, i
    integer :: size_world, rank_world

    integer, dimension(:), allocatable :: color_model_rec ! task_id
    integer, dimension(:), allocatable :: color_couple_rec
    integer, dimension(:), allocatable :: rank_model_rec
    integer, dimension(:), allocatable :: rank_filter_rec
    integer, dimension(:), allocatable :: rank_couple_rec
    integer, dimension(:), allocatable :: rank_world_rec

    logical, dimension(:), allocatable :: filterpe_rec

    integer, parameter    :: root = 0

    CALL MPI_Comm_size(MPI_COMM_WORLD, size_world, MPIerr)
    CALL MPI_Comm_rank(MPI_COMM_WORLD, rank_world, MPIerr)

    if( rank_world == root ) then
      allocate(color_model_rec(size_world))
      allocate(color_couple_rec(size_world))
      allocate(rank_model_rec(size_world))
      allocate(rank_filter_rec(size_world))
      allocate(rank_couple_rec(size_world))
      allocate(rank_world_rec(size_world))
      allocate(filterpe_rec(size_world))
    else
      allocate(color_model_rec(0))
      allocate(color_couple_rec(0))
      allocate(rank_model_rec(0))
      allocate(rank_filter_rec(0))
      allocate(rank_couple_rec(0))
      allocate(rank_world_rec(0))
      allocate(filterpe_rec(0))
    end if

    call MPI_Gather (task_id,     1, MPI_INTEGER, &
                    color_model_rec, 1, MPI_INTEGER, &
                    root,                   &
                    MPI_COMM_WORLD, MPIerr)

    call MPI_Gather (color_couple,    1, MPI_INTEGER, &
                    color_couple_rec, 1, MPI_INTEGER, &
                    root,                   &
                    MPI_COMM_WORLD, MPIerr)

    call MPI_Gather (mytid,          1, MPI_INTEGER, &
                    rank_model_rec, 1, MPI_INTEGER, &
                    root,                   &
                    MPI_COMM_WORLD, MPIerr)

    call MPI_Gather (rank_filter,     1, MPI_INTEGER, &
                    rank_filter_rec, 1, MPI_INTEGER, &
                    root,                   &
                    MPI_COMM_WORLD, MPIerr)

    call MPI_Gather (rank_couple,     1, MPI_INTEGER, &
                    rank_couple_rec, 1, MPI_INTEGER, &
                    root,                   &
                    MPI_COMM_WORLD, MPIerr)

    call MPI_Gather (rank_world,     1, MPI_INTEGER, &
                    rank_world_rec, 1, MPI_INTEGER, &
                    root,                   &
                    MPI_COMM_WORLD, MPIerr)

    call MPI_Gather (filterpe,     1, MPI_INTEGER, &
                    filterpe_rec, 1, MPI_INTEGER, &
                    root,                   &
                    MPI_COMM_WORLD, MPIerr)

    CALL MPI_Barrier(MPI_COMM_WORLD, MPIerr)
    if( rank_world == root ) then
      WRITE (*, '(/18x, a)') 'PE configuration:'
      WRITE (*, '(2x, a6, a9, a10, a14, a13, /2x, a5, a9, a7, a7, a7, a7, a7, /2x, a)') &
          'world', 'filter', 'model', 'couple', 'filterPE', &
          'rank', 'rank', 'task', 'rank', 'task', 'rank', 'T/F', &
          '----------------------------------------------------------'

      do i=1, size_world
          if( filterpe_rec(i) .eqv. .true.) then
            WRITE (*, '(2x, i4, 4x, i4, 4x, i3, 4x, i3, 4x, i3, 4x, i3, 5x, l3)') &
            rank_world_rec(i), rank_filter_rec(i), color_model_rec(i), rank_model_rec(i), &
            color_couple_rec(i), rank_couple_rec(i), filterpe_rec(i)
          else
            WRITE (*,'(2x, i4, 12x, i3, 4x, i3, 4x, i3, 4x, i3, 5x, l3)') &
            rank_world_rec(i), color_model_rec(i), rank_model_rec(i), &
            color_couple_rec(i), rank_couple_rec(i), filterpe_rec(i)
          end if
      end do

      WRITE (*, '(/a)') ''

    end if

    deallocate(color_model_rec)
    deallocate(rank_model_rec)
    deallocate(rank_filter_rec)
    deallocate(rank_couple_rec)
    deallocate(rank_world_rec)
    deallocate(filterpe_rec)

    CALL MPI_Barrier(MPI_COMM_WORLD, MPIerr)

  end subroutine print_mpi_layout

  !> Gathers the world ranks of all processes in this model run's communicator to its local rank 0 and prints the model-run-id to world-rank mapping.
  subroutine print_modelrun_world_rank_mapping

    implicit none

    integer :: size_world, rank_world
    integer :: size_model, rank_model
    integer :: MPIerr

    integer, dimension(:), allocatable :: rank_world_rec

    CALL MPI_Comm_size(MPI_COMM_WORLD, size_world, MPIerr)
    CALL MPI_Comm_rank(MPI_COMM_WORLD, rank_world, MPIerr)

    CALL MPI_Comm_size(TIEGCM_WORLD, size_model, MPIerr)
    CALL MPI_Comm_rank(TIEGCM_WORLD, rank_model, MPIerr)


    if(mytid == 0) then
        allocate(rank_world_rec(size_model))
    end if

    call MPI_Gather (rank_world,     1, MPI_INTEGER, &
                    rank_world_rec, 1, MPI_INTEGER, &
                    0,                   &
                    TIEGCM_WORLD, MPIerr)

    if(mytid == 0) then
        write(*, '(a, i4, a, 16i4)' ) 'model run id: ', task_id, '       associated world rank: ', rank_world_rec
        deallocate(rank_world_rec)
    end if

  end subroutine print_modelrun_world_rank_mapping

END MODULE print_parallel_info
