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
! PDAF pre/post-step hook: called before/after each analysis to inspect the ensemble and update calibration parameters.

!BOP
!
! !ROUTINE: prepoststep_ens_pdaf --- Used-defined Pre/Poststep routine for PDAF
!
! !INTERFACE:
SUBROUTINE prepoststep_ens_pdaf(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
! 
! The routine is called for global filters (e.g. SEIK)
! before the analysis and after the ensemble transformation.
! For local filters (e.g. LSEIK) the routine is called
! before and after the loop over all local analysis
! domains.
! The routine provides full access to the state 
! estimate and the state ensemble to the user.
! Thus, user-controlled pre- and poststep 
! operations can be performed here. For example 
! the forecast and the analysis states and ensemble
! covariance matrix can be analyzed, e.g. by 
! computing the estimated variances. 
! For the offline mode, this routine is the place
! in which the writing of the analysis ensemble
! can be performed.
!
! If a user considers to perform adjustments to the 
! estimates (e.g. for balances), this routine is 
! the right place for it.
!
! called by all filter processes
!
! !USES:

  ! pdaf
  use pdaf, only: PDAF_get_localfilter

  ! intern
  use configuration, only: cfg_calibration
  use mod_assimilation, only: co_estimate_dynamics_this_step
  use print_parallel_info, only: print_ranks

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step        ! Current time step (negative for call after forecast)
  INTEGER, INTENT(in) :: dim_p       ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens     ! Size of state ensemble
  INTEGER, INTENT(in) :: dim_ens_p   ! PE-local size of ensemble
  INTEGER, INTENT(in) :: dim_obs_p   ! PE-local dimension of observation vector
  REAL, INTENT(inout) :: state_p(dim_p) ! PE-local forecast/analysis state
  ! The array 'state_p' is not generally not initialized in the case of SEIK.
  ! It can be used freely here.
  REAL, INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1) ! Inverse of matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)      ! PE-local state ensemble
  INTEGER, INTENT(in) :: flag        ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: PDAF_get_state      (as U_prepoststep)
! Called by: PDAF_X_update       (as U_prepoststep)
! Calls: MPI_send
! Calls: MPI_recv
!EOP

! *** local variables ***
  LOGICAL, SAVE :: firsttime = .TRUE. ! Routine is called for first time?
  CHARACTER(len=3),dimension(3) :: anastr = (/'ini','for','ana'/)  ! String for call type (initial, forecast, analysis)
  integer :: callType = -1

  REAL, ALLOCATABLE, dimension(:,:), save :: ens_forecast_p

!   REAL, ALLOCATABLE :: variance_p(:)   ! local model state variances
!   REAL, ALLOCATABLE :: mean_p(:)       ! local model state mean

  INTEGER :: localfilter          ! Flag for domain-localized filter (1=true)


! **********************
! *** INITIALIZATION ***
! **********************


  IF (firsttime) THEN
    WRITE (*, '(8x, a)') 'Analyze initial state ensemble'
    callType = 1
  ELSE
    IF (step<0) THEN
        WRITE (*, '(8x, a)') 'Analyze and write forecasted state ensemble'
        callType = 2
    ELSE
        WRITE (*, '(8x, a)') 'Analyze and write assimilated state ensemble'
        callType = 3
    END IF
  END IF

  WRITE(*,*) 'prepost ', anastr(callType), ' step:', step, dim_p, dim_ens, dim_ens_p, dim_obs_p


! **********************************
! *** ensemble mean and variance ***
! **********************************

!   ALLOCATE(variance_p(dim_p))
!   ALLOCATE(mean_p(dim_p))
!   call compute_mean_and_variance_vec(ens_p, mean=mean_p, variance=variance_p, dim=2)

  CALL PDAF_get_localfilter(localfilter)

  if(localfilter==1) then
    if(cfg_calibration%apply) then
      ! **********************************
      ! *** store forecast ensemble    ***
      ! **********************************
      if(callType == 2)then
        if( allocated(ens_forecast_p)) deallocate(ens_forecast_p)
        allocate(ens_forecast_p, source=ens_p)
      end if

      ! **********************
      ! calibration parameters
      ! **********************
      if(callType == 3) then
        if(co_estimate_dynamics_this_step) then
          call update_calibration_parameters(dim_p, dim_ens, ens_p, ens_forecast_p)
        end if
      end if
    end if
  end if

! ********************
! *** finishing up ***
! ********************

  firsttime = .FALSE.

!    IF (ALLOCATED(variance_p)) DEALLOCATE(variance_p)
!    IF (ALLOCATED(mean_p)) DEALLOCATE(mean_p)

  ! Deallocate observation arrays
  CALL deallocate_obs_pdafomi(step)

END SUBROUTINE prepoststep_ens_pdaf


subroutine update_calibration_parameters(dim_p,dim_ens,ens_analysis_p, ens_forecast_p)

  ! tie-gcm
  use mpi_module, only: mytid

  ! intern
  use ensemble_module, only: update_calibration_parameters_lin_reg
  use state_module,  only: state_vector

  implicit none

  ! arguments
  integer, intent(in) :: dim_p, dim_ens
  real, dimension(dim_p, dim_ens), intent(inout) :: ens_analysis_p ! local augumented ensemble at analysis
  real, dimension(dim_p, dim_ens), intent(inout) :: ens_forecast_p ! local augumented ensemble at forecast

  ! local
  integer :: idx_state_0, idx_state_1, dim_state_p
  integer :: idx_cal_0, idx_cal_1, dim_cal_p


  idx_state_0 = state_vector%map%idx_R(state_vector%idx_f3d_0,mytid)%begin_p
  idx_state_1 = state_vector%map%idx_R(state_vector%idx_f3d_1,mytid)%back_p

  dim_state_p = (idx_state_1-idx_state_0)+1

  idx_cal_0 = state_vector%map%idx_R(state_vector%idx_cal_0,mytid)%begin_p
  idx_cal_1 = state_vector%map%idx_R(state_vector%idx_cal_1,mytid)%back_p

  dim_cal_p = (idx_cal_1-idx_cal_0)+1

  call update_calibration_parameters_lin_reg(dim_state_p, dim_cal_p, dim_ens, ens_analysis_p, ens_forecast_p)

end subroutine
