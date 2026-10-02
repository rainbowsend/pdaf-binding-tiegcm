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
! Calls PDAF's assimilation (analysis) step once per forecast cycle, dispatching to the local or global filter variant.

!>  Routine to call PDAF for analysis step
!!
!! This routine is called during the model integrations at each time 
!! step. It calls the filter-speific assimilation routine of PDAF 
!! (PDAF_assimilate_X), which checks whether the forecast phase is
!! completed. If so, the analysis step is computed inside PDAF
!!
SUBROUTINE assimilate_pdaf()

  ! pdaf
  use PDAF, &   ! Interface definitions to PDAF core routines
      only: PDAF3_assimilate_local_nondiagR, &
      PDAF3_assimilate_global, &
      PDAFomi_assimilate_lenkf, PDAFomi_generate_obs, PDAF_get_localfilter
  use PDAF_mod_core, &
      only: cnt_steps, nsteps

  ! intern
  use mod_parallel_pdaf, &        ! Parallelization variables
      only: rank_world, abort_parallel
!   use mod_assimilation, &         ! Filter variables
!       only: filtertype
  use configuration,&
      only: cfg_filter, cfg_calibration
  use mod_assimilation,&
      only: analysis_step_count, analysis_step_dynamics_count, &
            co_estimate_dynamics_this_step


  IMPLICIT NONE

! *** Local variables ***
  INTEGER :: status_pdaf          ! PDAF status flag
  INTEGER :: localfilter          ! Flag for domain-localized filter (1=true)

! External subroutines
!   (subroutine names are passed over to PDAF in the calls to 
!   PDAF_get_state and PDAF_assimilate_X. This allows the user 
!   to specify the actual name of a routine. However, the 
!   PDAF-internal name of a subroutine might be different from
!   the external name!)

  ! Interface between model and PDAF, and prepoststep
  EXTERNAL :: collect_state_pdaf, &   ! Collect a state vector from model fields
       distribute_state_pdaf, &       ! Distribute a state vector to model fields
       next_observation_pdaf, &       ! Provide time step of next observation
       prepoststep_ens_pdaf           ! User supplied pre/poststep routine
  ! Localization of state vector
  EXTERNAL :: init_n_domains_pdaf, &  ! Provide number of local analysis domains
       init_dim_l_pdaf                ! Initialize state dimension for local analysis domain
  ! Interface to PDAF-OMI for local and global filters
  EXTERNAL :: init_dim_obs_pdafomi, & ! Get dimension of full obs. vector for PE-local domain
       obs_op_pdafomi, &              ! Obs. operator for full obs. vector for PE-local domain
       init_dim_obs_l_pdafomi, &      ! Get dimension of obs. vector for local analysis domain
       localize_covar_pdafomi, &      ! Apply localization to covariance matrix in LEnKF
       prodRinvA_l_pdafomi
! ! Subroutine used for generating observations
!   EXTERNAL :: get_obs_f_pdaf          ! Get vector of synthetic observations from PDAF

  ! check whether at this time step the anlysis step is performed
  co_estimate_dynamics_this_step = .false.
  IF (cnt_steps+1 == nsteps) THEN
    analysis_step_count = analysis_step_count + 1
    if (cfg_calibration%apply) then
      if(MODULO(analysis_step_count, cfg_calibration%every)==0)then
        co_estimate_dynamics_this_step = .true.
        analysis_step_dynamics_count = analysis_step_dynamics_count + 1
      end if
    end if
  end if

  if(cfg_filter%open_loop .eqv. .false.) then

! *********************************
! *** Call assimilation routine ***
! *********************************

  ! Check  whether the filter is domain-localized
  CALL PDAF_get_localfilter(localfilter)

  ! Call assimilate routine for global or local filter
  IF (localfilter==1) THEN
        CALL PDAF3_assimilate_local_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdafomi, prodRinvA_l_pdafomi, &
          prepoststep_ens_pdaf, next_observation_pdaf, status_pdaf)
  ELSE
!      IF (filtertype==8) THEN
!         ! LEnKF has its own OMI interface routine
!         CALL PDAFomi_assimilate_lenkf(collect_state_pdaf, distribute_state_pdaf, &
!              init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_ens_pdaf, &
!              localize_covar_pdafomi, next_observation_pdaf, status_pdaf)
!      ELSE IF (filtertype==11) THEN
!         ! Observation generation has its own OMI interface routine
!         CALL PDAFomi_generate_obs(collect_state_pdaf, distribute_state_pdaf, &
!              init_dim_obs_pdafomi, obs_op_pdafomi, get_obs_f_pdaf, &
!              prepoststep_ens_pdaf, next_observation_pdaf, status_pdaf)
!      ELSE
!         ! All global filters except LEnKF
        CALL PDAF3_assimilate_global(collect_state_pdaf, distribute_state_pdaf, &
             init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_ens_pdaf, &
             next_observation_pdaf, status_pdaf)
!      END IF
  END IF

  ! Check for errors during execution of PDAF

  IF (status_pdaf /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in PDAF_put_state - stopping! (PE ', rank_world,')'
     CALL  abort_parallel()
  END IF
  end if

END SUBROUTINE assimilate_pdaf
