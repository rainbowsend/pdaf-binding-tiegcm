! Based on PDAF template/tutorial code.
! Copyright (c) 2004-2026 Lars Nerger, Alfred Wegener Institute,
! Helmholtz Center for Polar and Marine Research, Bremerhaven, Germany.
!
! Modified for the TIE-GCM/PDAF coupling.
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
! Module of shared MPI parallelization variables/communicators used across PDAF routines.

!BOP
!
! !MODULE:
MODULE mod_parallel_pdaf

! !DESCRIPTION:
! This modules provides variables for the MPI parallelization
! to be shared between model-related routines. The are variables
! that are used in the model, even without PDAF and additional
! variables that are only used, if data assimialtion with PDAF
! is performed.
! In addition methods to initialize and finalize MPI are provided.
! The initialization routine is only for the model itself, the 
! more complex initialization of communicators for xecution with
! PDAF is peformed in init\_parallel\_pdaf.
!
! !USES:

  ! extern
  use mpi_f08

  IMPLICIT NONE
  SAVE 

! !PUBLIC DATA MEMBERS:
  ! Additional variables for use with PDAF
  INTEGER :: n_modeltasks ! Number of parallel model tasks ! set in config file
  INTEGER :: task_id     ! Index of my model task (1,...,n_modeltasks) color of MPI_COMM
  INTEGER :: local_ntask ! number of processing elements per ensemble
  TYPE(MPI_Comm) :: COMM_filter ! MPI communicator for filter PEs
  INTEGER :: rank_filter
  INTEGER :: npes_filter ! number of processing elements in communicator COMM_filter
  LOGICAL :: filterpe    ! Whether we are on a PE in a COMM_filter
  TYPE(MPI_Comm) :: COMM_couple ! MPI communicator for coupling filter and model
  INTEGER :: color_couple
  INTEGER :: rank_couple
  INTEGER :: npes_couple ! number of processing elements in communicator COMM_couple
  INTEGER :: rank_world
  INTEGER :: npes_world  ! number of processing elements in communicator MPI_COMM_world
  INTEGER :: MPIerr      ! Error flag for MPI
  TYPE(MPI_Status) :: MPIstatus ! Status array for MPI
!EOP

  contains

  subroutine abort_parallel
    IMPLICIT NONE
    integer :: ier
    write(*,'(a)') 'PDAF is aborting MPI'
    call MPI_Abort(MPI_COMM_WORLD, 1, ier)
    ! shutdown function
  end subroutine abort_parallel

END MODULE mod_parallel_pdaf
