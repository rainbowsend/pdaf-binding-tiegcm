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
! Writes the PDAF state vector back into TIE-GCM model fields and applies constraints and parameter updates.

!BOP
!
! !ROUTINE: distribute_state_pdaf --- Initialize model fields from state vector
!
! !INTERFACE:
SUBROUTINE distribute_state_pdaf(dim_p, state_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! During the forecast phase of the filter this
! subroutine is called from PDAF\_get\_state
! supplying a model state which has to be evolved. 
! The routine has to initialize the fields of the 
! model (typically available through a module) from 
! the state vector of PDAF. With parallelization, 
! MPI communication might be required to 
! initialize all subdomains on the model PEs.
!
! The routine is executed by each process that is
! participating in the model integrations.
!
! For the dummy model and PDAF with domain
! decomposition the state vector and the model
! field are identical. Hence, the field array 
! is directly initialized from an ensemble 
! state vector by each model PE.
!

  ! extern
  use mpi_f08

  ! tie-gcm
  use fields_module,&
      only: f4d, itc
  use mpi_module,&
      only: mytid, handle_mpi_err, mp_periodic_f4d

  ! intern
  use array_print_module,&
      only: printMat
  use configuration,&
      only: cfg_calibration, cfg_output
 use model_parameter_IO_module,&
      only: model_parameter_writer
  use constraints,&
      only: constrain_state
  use mod_assimilation,&
      only: analysis_step_count, co_estimate_dynamics_this_step
  use result_writer_frontend,&
      only: compute_and_write_results
  use state_module,&
      only: state_vector, levX0, &
            levX1, nlevX, idx_intern


!! For test output
!   use mod_parallel_pdaf, only: task_id
!   use state_module, only:  fd_name
!   use mpi_module, only:ntask

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p           ! PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  ! PE-local state vector

! !CALLING SEQUENCE:
! Called by: PDAF_get_state      (as U_dist_state)
! Called by: PDAF_assimilate_X   (as U_coll_state)
!EOP

! *** local variables ***
  INTEGER :: i        ! Counters

! For test output
!   INTEGER :: ier
!   INTEGER :: j

! *******************************************
! *** Initialize model fields from state  ***
!********************************************

    write(*,*) "Initialize model from state (", analysis_step_count,")"

    do i= state_vector%idx_f3d_0, state_vector%idx_f3d_1
        f4d(state_vector%fd_idx(i))%data(           levX0 : levX1,                    &
                                   idx_intern(mytid)%lon0 : idx_intern(mytid)%lon1,   &
                                   idx_intern(mytid)%lat0 : idx_intern(mytid)%lat1,   &
                                                          itc                  )      &
        = reshape( &
           state_p( state_vector%map%idx_R(i, mytid)%begin_p :&
                    state_vector%map%idx_R(i, mytid)%back_p), &
           (/nlevX, idx_intern(mytid)%nlons, idx_intern(mytid)%nlats /) &
          )
    end do


    ! 'mp_periodic_f4d' is called in advance.F after 'assimilate_pdaf' anyway, no need to call it here
    ! call mp_periodic_f4d(itc)

    if((analysis_step_count>0) .and. (cfg_output%save_unconstrained_analysis)) then
      call compute_and_write_results("unconstrained_analysis")
    end if

    ! PDAF_get_state in init_pdaf calls this subroutine.
    ! At that time (analysis_step_count==0) no analysis step was perfomred before
    if(analysis_step_count>0)then
      call constrain_state()
    end if

    ! barm is prognostic field (required to initialize TIE-GCM), it is calculated in "addiag" which called in "advance" routine. Before that computation barm is not used -> we do not need to update it


! *******************************************
! *** init model parameters from state    ***
!********************************************

    ! must be called at analysis_step_count==0 to initalize the ensemble and to write initail dynamics
    if(cfg_calibration%apply .and. (co_estimate_dynamics_this_step .or. analysis_step_count==0)) then
      call state_vector%distribute_calibration(state_p)
      call model_parameter_writer%write(state_p)
    end if

! *******************************************
! *** write analysis step to file         ***
!********************************************

    ! PDAF_get_state in init_pdaf calls this subroutine.
    ! At that time (analysis_step_count==0) no analysis step was perfomred before
    if(analysis_step_count>0)then
      call compute_and_write_results("analysis")
    end if

!     ! write fields to console to check MPI scatter and reshape
!     call MPI_Barrier(MPI_COMM_WORLD,ier)
!
!     if( task_id == 1 ) then
!         do j =0, ntask-1
!             call MPI_Barrier(COMM_model,ier)
!             if(mytid == j) then
!                 do i =1, nfields
!                     call printMat(f4d(fd_idx(i))%data(:,:,:,itc), 'state_p_field ' // fd_name(i))
!                 end do
!             end if
!             call MPI_Barrier(COMM_model,ier)
!         end do
!     end if
!
!     call MPI_Barrier(MPI_COMM_WORLD,ier)

END SUBROUTINE distribute_state_pdaf
