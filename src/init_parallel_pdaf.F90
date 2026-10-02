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
! Sets up MPI communicators (filter/model/coupling) for running PDAF alongside TIE-GCM.

!BOP
!
! !ROUTINE: init_parallel_pdaf --- Initialize communicators for PDAF
!
! !INTERFACE:
SUBROUTINE init_parallel_pdaf(dim_ens, screen)

! !DESCRIPTION:
! Parallelization routine for a model with 
! attached PDAF. The subroutine is called in 
! the main program subsequently to the 
! initialization of MPI. It initializes
! MPI communicators for the model tasks, filter 
! tasks and the coupling between model and
! filter tasks. In addition some other variables 
! for the parallelization are initialized.
! The communicators and variables are handed
! over to PDAF in the call to 
! PDAF\_filter\_init.
!
! 3 Communicators are generated:\\
! - COMM\_filter: Communicator in which the
!   filter itself operates\\
! - COMM\_model: Communicators for parallel
!   model forecasts\\
! - COMM\_couple: Communicator for coupling
!   between models and filter\\
! Other variables that have to be initialized are:\\
! - filterpe - Logical: Does the PE execute the 
! filter?\\
! - my\_ensemble - Integer: The index of the PE's 
! model task\\
! - local\_npes\_model - Integer array holding 
! numbers of PEs per model task
!
! For COMM\_filter and COMM\_model also
! the size of the communicators (npes\_filter and 
! npes\_model) and the rank of each PE 
! (mype\_filter, mype\_model) are initialized. 
! These variables can be used in the model part 
! of the program, but are not handed over to PDAF.
!
! This variant is for a domain decomposed 
! model.
!
! NOTE: 
! This is a template that is expected to work 
! with many domain-decomposed models. However, 
! it might be necessary to adapt the routine 
! for a particular model. Inportant is that the
! communicator COMM_model equals the communicator
! used in the model. If one plans to run a parallel 
! ensemble forecast (that is using multiple model
! tasks), COMM_model cannot be MPI_COMM_WORLD! Thus,
! if the model uses MPI_COMM_WORLD it has to be
! replaced by an alternative communicator named,
! e.g., COMM_model.
!
! !USES:

  ! extern
  use mpi_f08

  ! tie-gcm
  use mpi_module, only: TIEGCM_WORLD, handle_mpi_err, ntask, mytid, tasks

  ! intern
  use mod_parallel_pdaf, &
         only: rank_filter, npes_filter, COMM_filter, filterpe, &
         n_modeltasks, local_ntask, task_id, COMM_couple, &
         MPI_COMM_WORLD, rank_world, npes_world, MPIerr,  &
         rank_couple, npes_couple, color_couple
  USE print_parallel_info, &
        only: print_modelrun_world_rank_mapping, print_mpi_layout

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(inout) :: dim_ens ! Ensemble size or number of EOFs (only SEEK)
  ! Often dim_ens=0 when calling this routine, because the real ensemble size
  ! is initialized later in the program. For dim_ens=0 no consistency check
  ! for ensemble size with number of model tasks is performed.
  INTEGER, INTENT(in)    :: screen ! Whether screen information is shown

! !CALLING SEQUENCE:
! Called by: tgcm.F
! Calls: MPI_Comm_size
! Calls: MPI_Comm_rank
! Calls: MPI_Comm_split
! Calls: MPI_Barrier
!EOP

  ! local variables
  INTEGER :: i, j                 ! Counters
  TYPE(mpi_comm) :: COMM_ensemble ! Communicator of all PEs doing model tasks
  INTEGER :: mype_ens, npes_ens   ! rank and size in COMM_ensemble
  INTEGER :: pe_index             ! Index of PE
  INTEGER :: my_color             ! Variables for communicator-splitting
  LOGICAL :: iniflag              ! Flag whether MPI is initialized

  ! copy from mpi.F mp_init
  integer,dimension(12) :: valid_pecounts=(/1,4,6,8,12,16,24,32,48,64,72,80/)


  ! *** Initialize MPI if not yet initialized ***
  CALL MPI_Initialized(iniflag, MPIerr)
  IF (.not.iniflag) THEN
     CALL MPI_Init(MPIerr)
     if (MPIerr /= MPI_SUCCESS) call  handle_mpi_err(MPIerr,'init pdaf mpi')
  END IF

  ! *** Initialize PE information on COMM_world ***
  ! save number of tasks and rank of MPI_COMM_WORLD to 
  ! variables in mod_parallel_pdaf
  CALL MPI_Comm_size(MPI_COMM_WORLD, npes_world, MPIerr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank_world, MPIerr)


  ! *** Initialize communicators for ensemble evaluations ***
  IF (rank_world == 0) &
       WRITE (*, '(/1x, a)') 'Initialize communicators for assimilation with PDAF'


  ! *** Check consistency of number of parallel ensemble tasks ***
  consist1: IF (n_modeltasks > npes_world) THEN
     ! *** # parallel tasks is set larger than available PEs ***
     n_modeltasks = npes_world
     IF (rank_world == 0) WRITE (*, '(3x, a)') &
          '!!! Resetting number of parallel ensemble tasks to total number of PEs!'
  END IF consist1
  IF (dim_ens > 0) THEN
     ! Check consistency with ensemble size
     consist2: IF (n_modeltasks > dim_ens) THEN
        ! # parallel ensemble tasks is set larger than ensemble size
        n_modeltasks = dim_ens
        IF (rank_world == 0) WRITE (*, '(5x, a)') &
             '!!! Resetting number of parallel ensemble tasks to number of ensemble states!'
     END IF consist2
  END IF


  ! ***              COMM_ENSEMBLE                ***
  ! *** Generate communicator for ensemble runs   ***
  ! *** only used to generate model communicators ***
  call MPI_Comm_dup(MPI_COMM_WORLD, COMM_ensemble, MPIerr)
  if (MPIerr /= MPI_SUCCESS) call  handle_mpi_err(MPIerr,'init pdaf dublicate world')

  CALL MPI_Comm_Size(COMM_ensemble, npes_ens, MPIerr)
  CALL MPI_Comm_Rank(COMM_ensemble, mype_ens, MPIerr)

  ! *** Store # PEs per ensemble                 ***
  ! *** used for info on PE 0 and for generation ***
  ! *** of model communicators on other Pes      ***


  local_ntask = FLOOR(REAL(npes_world) / REAL(n_modeltasks))

  if(  .not.any(local_ntask == valid_pecounts) ) then
    call shutdown( '(number of processes / ensemble size) must be in [1,4,6,8,12,16,24,32,48,64,72,80]. See mpi.F, mp_init ' )
  end if

  IF (rank_world == 0) &
    write (*,*) 'n global task ', npes_world, ' n_modeltasks', n_modeltasks, ' n local task', local_ntask


  ! ***              COMM_MODEL               ***
  ! *** Generate communicators for model runs ***
  ! *** (Split COMM_ENSEMBLE)                 ***
  pe_index = 0
  doens1: DO i = 1, n_modeltasks
     DO j = 1, local_ntask
        IF (mype_ens == pe_index) THEN
           task_id = i
           EXIT doens1
        END IF
        pe_index = pe_index + 1
     END DO
  END DO doens1

  CALL MPI_Comm_split(COMM_ensemble, task_id, mype_ens, &
       TIEGCM_WORLD, MPIerr)
  if (MPIerr /= MPI_SUCCESS) call  handle_mpi_err(MPIerr,'split model')

  ! *** Re-initialize PE informations   ***
  ! *** according to model communicator ***
  CALL MPI_Comm_Size(TIEGCM_WORLD, ntask, MPIerr)
  CALL MPI_Comm_Rank(TIEGCM_WORLD, mytid, MPIerr)

  if (screen > 1) then
    write (*,*) 'MODEL: mype(w)= ', rank_world, '; model task: ', task_id, &
         '; mype(m)= ', mytid, '; npes(m)= ', ntask
  end if


  ! Init flag FILTERPE (all PEs of model task 1)
  IF (task_id == 1) THEN
     filterpe = .TRUE.
  ELSE
     filterpe = .FALSE.
  END IF

  ! ***         COMM_FILTER                 ***
  ! *** Generate communicator for filter    ***
  ! *** For simplicity equal to COMM_couple ***
  my_color = task_id

  CALL MPI_Comm_split(MPI_COMM_WORLD, my_color, rank_world, &
       COMM_filter, MPIerr)
  if (MPIerr /= MPI_SUCCESS) call  handle_mpi_err(MPIerr,'split filter')

  ! *** Initialize PE informations         ***
  ! *** according to coupling communicator ***
  CALL MPI_Comm_Size(COMM_filter, npes_filter, MPIerr)
  CALL MPI_Comm_Rank(COMM_filter, rank_filter, MPIerr)


  ! ***              COMM_COUPLE                 ***
  ! *** Generate communicators for communication ***
  ! *** between model and filter PEs             ***
  ! *** (Split COMM_ENSEMBLE)                    ***

  color_couple = rank_filter + 1

  CALL MPI_Comm_split(MPI_COMM_WORLD, color_couple, rank_world, &
       COMM_couple, MPIerr)
  if (MPIerr /= MPI_SUCCESS) call  handle_mpi_err(MPIerr,'split couple')

  ! *** Initialize PE informations         ***
  ! *** according to coupling communicator ***
  CALL MPI_Comm_Size(COMM_couple, npes_couple, MPIerr)
  CALL MPI_Comm_Rank(COMM_couple, rank_couple, MPIerr)


  IF (screen > 0) THEN
    ! write ordered (only world root writes) summary 
    call print_mpi_layout
    call print_modelrun_world_rank_mapping
  end if

! ******************************************************************************
! *** Initialize model equivalents to COMM_model, ntask, and mytid ***
! ******************************************************************************

  ! If the names of the variables for COMM_model, ntask, and
  ! mytid are different in the numerical model, the 
  ! model-internal variables should be initialized at this point.

  ! Template reminder - delete when implementing functionality
  !WRITE (*,*) 'TEMPLATE init_parallel_pdaf.F90: Initialize model communicator here!'

  if (allocated (tasks)) deallocate (tasks)
  allocate(tasks(0:ntask-1))

END SUBROUTINE init_parallel_pdaf
