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
! Dispatches PDAF-OMI observation call-backs (dimensions, obs. operator, R^-1 product) to all registered observation types.

!> callback_obs_pdafomi
!!
!! This file provides interface routines between the call-back routines
!! of PDAF and the observation-specific routines in PDAF-OMI. This structure
!! collects all calls to observation-specific routines in this single file
!! to make it easier to find the routines that need to be adapted.
!!
!! The routines here are mainly pure pass-through routines. Thus they
!! simply call one of the routines from PDAF-OMI. Partly some addtional
!! variable is required, e.g. to specify the offset of an observation
!! in the observation vector containing all observation types. These
!! cases are described in the routines.
!!
!! **Adding an observation type:**
!!   When adding an observation type, one has to add one module
!!   obs_TYPE_pdafomi (based on the template obs_TYPE_pdafomi_TEMPLATE.F90).
!!   In addition one has to add a call to the different routines include
!!   in this file. It is recommended to keep the order of the calls
!!   consistent over all files. 
!! 
!-------------------------------------------------------------------------------

!> Call-back routine for init_dim_obs
!!
!! This routine calls the observation-specific
!! routines init_dim_obs_TYPE.
!!
SUBROUTINE init_dim_obs_pdafomi(step, dim_obs)

  ! Include functions for different observations

  ! intern
  use observations_module, only: observations_assim

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step     !< Current time step
  INTEGER, INTENT(out) :: dim_obs  !< Dimension of full observation vector

! *** Local variables ***

  integer :: dim_obs_tmp
  integer :: i


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

  ! Initialize number of observations
  dim_obs = 0

  ! Call observation-specific routines
  ! The routines are independent, so it is not relevant
  ! in which order they are called

  do i=1,size(observations_assim)
    if(observations_assim(i)%ptr%is_available_at_next_assim_step) then
      call observations_assim(i)%ptr%init_dim_obs(step, dim_obs_tmp)
      dim_obs = dim_obs + dim_obs_tmp
    end if
  end do

END SUBROUTINE init_dim_obs_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for obs_op
!!
!! This routine calls the observation-specific
!! routines obs_op_TYPE.
!!
SUBROUTINE obs_op_pdafomi(step, dim_p, dim_obs, state_p, ostate)

  ! intern
  use observations_module, only: observations_assim

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step                 !< Current time step
  INTEGER, INTENT(in) :: dim_p                !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_obs              !< Dimension of full observed state
  REAL, INTENT(in)    :: state_p(dim_p)       !< PE-local model state
  REAL, INTENT(inout) :: ostate(dim_obs)      !< PE-local full observed state

  integer :: i


! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

  ! The order of these calls is not relevant as the setup
  ! of the overall observation vector is defined by the
  ! order of the calls in init_dim_obs_pdafomi

  do i=1,size(observations_assim)
    if(observations_assim(i)%ptr%is_available_at_next_assim_step) then
      call observations_assim(i)%ptr%obs_op(dim_p, dim_obs, state_p, ostate)
    end if
  end do

END SUBROUTINE obs_op_pdafomi

!-------------------------------------------------------------------------------
!> Call-back routine for prodRinvA_l
!!
!! This routine calls the observation-specific
!! routines prodRinvA_l.
!!
SUBROUTINE prodRinvA_l_pdafomi(domain_p, step, dim_obs_l, rank, obs_l, A_l, C_l)

  ! intern
  use observations_module, only: observations_assim

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: domain_p          !< Index of current local analysis domain
  INTEGER, INTENT(in) :: step              !< Current time step
  INTEGER, INTENT(in) :: dim_obs_l         !< Dimension of local observation vector
  INTEGER, INTENT(in) :: rank              !< Rank of initial covariance matrix
  REAL, INTENT(in)    :: obs_l(dim_obs_l)  !< Local vector of observations
  REAL, INTENT(inout) :: A_l(dim_obs_l, rank) !< Input matrix
  REAL, INTENT(out)   :: C_l(dim_obs_l, rank) !< Output matrix

  ! local
  integer :: i

  ! *** Compute
  ! ***                  -1
  ! ***           C = W R   A

  ! The order of these calls is not relevant as the setup
  ! of the overall observation vector is defined by the
  ! order of the calls in init_dim_obs_pdafomi

  do i=1,size(observations_assim)
    if(observations_assim(i)%ptr%is_available_at_next_assim_step) then
      call observations_assim(i)%ptr%prodRinvA_l(domain_p, step, dim_obs_l, rank, obs_l, A_l, C_l)
    end if
  end do

END SUBROUTINE prodRinvA_l_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for init_dim_obs_l
!!
!! This routine calls the routine PDAFomi_init_dim_obs_l
!! for each observation type
!!
SUBROUTINE init_dim_obs_l_pdafomi(domain_p, step, dim_obs, dim_obs_l)

  ! intern
  use observations_module, only: observations_assim

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: domain_p   !< Index of current local analysis domain
  INTEGER, INTENT(in)  :: step       !< Current time step
  INTEGER, INTENT(in)  :: dim_obs    !< Full dimension of observation vector
  INTEGER, INTENT(out) :: dim_obs_l  !< Local dimension of observation vector

  integer :: i

! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

  ! Call init_dim_obs_l specific for each observation

  do i=1,size(observations_assim)
    if(observations_assim(i)%ptr%is_available_at_next_assim_step) then
      call observations_assim(i)%ptr%init_dim_obs_l(domain_p, step, dim_obs, dim_obs_l)
    end if
  end do

END SUBROUTINE init_dim_obs_l_pdafomi



! !-------------------------------------------------------------------------------
! !> Call-back routine for localize_covar
! !!
! !! This routine calls the routine PDAFomi_localize_covar
! !! for each observation type to apply covariance
! !! localization in the LEnKF.
! !!
! SUBROUTINE localize_covar_pdafomi(dim_p, dim_obs, HP_p, HPH)
! 
!   ! Include functions for different observations
!   USE obs_TYPE_pdafomi, ONLY: localize_covar_TYPE
! 
!   IMPLICIT NONE
! 
! ! *** Arguments ***
!   INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
!   INTEGER, INTENT(in) :: dim_obs               !< number of observations
!   REAL, INTENT(inout) :: HP_p(dim_obs, dim_p)  !< PE local part of matrix HP
!   REAL, INTENT(inout) :: HPH(dim_obs, dim_obs) !< Matrix HPH
! 
! ! *** local variables ***
!   REAL, ALLOCATABLE :: coords_p(:,:) ! Coordinates of PE-local state vector entries
! 
! 
! ! **********************
! ! *** INITIALIZATION ***
! ! **********************
! 
!   ! Initialize coordinate array
! 
!   ! One needs to provide the array COORDS_P holding the coordinates of each
!   ! element of the process-local state vector. Each column of the array holds
!   ! the information for one element. The array can be initialized here using
!   ! information on the model grid.
! 
!   ! ALLOCATE(coords_p(NROWS, dim_p))
! 
!   ! coords_p = ...
! 
! 
! ! *************************************
! ! *** Apply covariance localization ***
! ! *************************************
! 
!   ! Call localize_covar specific for each observation
!   CALL localize_covar_TYPE(dim_p, dim_obs, HP_p, HPH, coords_p)
! 
! 
! ! ****************
! ! *** Clean up ***
! ! ****************
! 
!   ! DEALLOCATE(coords_p)
! 
! END SUBROUTINE localize_covar_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for deallocate_obs
!!
!! This routine calls the routine PDAFomi_deallocate_obs
!! for each observation type
!!
SUBROUTINE deallocate_obs_pdafomi(step)

  ! pdaf
  USE PDAF, ONLY: PDAFomi_deallocate_obs

  ! intern
  use observations_module, only: observations_assim

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step   !< Current time step

  integer :: i

! *************************************
! *** Deallocate observation arrays ***
! *************************************
  write(*,*) 'Deallocate observation arrays'

  ! We need one call for each observation type

  do i=1,size(observations_assim)
    if(observations_assim(i)%ptr%is_available_at_next_assim_step) then
      call PDAFomi_deallocate_obs(observations_assim(i)%ptr%full_obs)
    end if
  end do

END SUBROUTINE deallocate_obs_pdafomi
