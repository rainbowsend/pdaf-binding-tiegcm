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
! PDAF callback that determines the number of model steps until the next available observation and whether assimilation should stop.

!BOP
!
! !ROUTINE: next_observation_pdaf --- Initialize information on next observation
!
! !INTERFACE:
SUBROUTINE next_observation_pdaf(stepnow, nsteps, doexit, time)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! The subroutine is called before each forecast phase
! by PDAF\_get\_state. It has to initialize the number 
! of time steps until the next available observation 
! (nsteps) and the current model time (time). In 
! addition the exit flag (exit) has to be initialized.
! It indicates if the data assimilation process is 
! completed such that the ensemble loop in the model 
! routine can be exited.
!
! The routine is called by all processes. 
!         
! Version for the dummy model. Identical for
! mode- and domain-decomposition .
!
! !USES:

  ! extern
  use esmf

  ! tie-gcm
  use input_module, only: step  ! TIE-GCM temporal resolution

  ! intern
  use configuration, only: cfg_filter
  use mod_parallel_pdaf, only: rank_world
  use observations_module, only: observations_assim
  use time_module, only: get_current_modeltime, timeinterval_to_tgcm_steps, model_end

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: stepnow  ! Number of the current time step
  INTEGER, INTENT(out) :: nsteps   ! Number of time steps until next obs
  INTEGER, INTENT(out) :: doexit   ! Whether to exit forecasting (1 for exit)
  REAL, INTENT(out)    :: time     ! Current model (physical) time

  ! local
  type(ESMF_Time) :: current_model_time, future_model_time
  type(ESMF_TimeInterval) :: forecast_duration
  type(ESMF_TimeInterval) :: time_till_first_analysis

  integer :: i

  character(len=80) :: current_timestr, future_timestr

  logical, save :: first_call = .true.

  logical, allocatable, dimension(:) :: obs_available
! !CALLING SEQUENCE:
! Called by: PDAF_get_state   (as U_next_obs)
!EOP


! *******************************************************
! *** Set number of time steps until next observation ***
! *******************************************************

  time = 0.0          ! Not used in this implementation

  ! --- current model time
  call get_current_modeltime(current_model_time)
  call ESMF_TimeGet(current_model_time, timeString=current_timestr)


  ! time the next observation is scheduled
  ! ( assume constant time step interval) 
  ! check whether observation exists
  ! there must be an observation at the exact time
  ! --> Avoid temporal interpolation

  call ESMF_TimeIntervalSet( forecast_duration, s = cfg_filter%forecast_duration_sec )

  if(first_call)then
    call ESMF_TimeIntervalSet( time_till_first_analysis, s = cfg_filter%first_analysis_step_sec )
    future_model_time = current_model_time + time_till_first_analysis
    first_call = .false.
  else
    future_model_time = current_model_time + forecast_duration
  end if

  allocate(obs_available(size(observations_assim)))

  search_next_obs: do while(future_model_time <= model_end)
    call ESMF_TimeGet(future_model_time, timeString=future_timestr)
    write(*,*) 'checking available observations at epoch '//trim(future_timestr)


    do i=1,size(observations_assim)
        obs_available(i) = observations_assim(i)%ptr%next_observation(future_model_time)
        observations_assim(i)%ptr%is_available_at_next_assim_step = obs_available(i)
    end do

    if (count(obs_available)>0) then
      exit search_next_obs
    end if

    future_model_time = future_model_time + forecast_duration
  end do search_next_obs

  if (count(obs_available)>0) then
    forecast_duration = future_model_time-current_model_time
    nsteps = timeinterval_to_tgcm_steps(forecast_duration)
    doexit = 0          ! Not used in this implementation

    IF (rank_world == 0) then
          WRITE (*, '(i7, 3a, 3x, a, i7, 3a)') &
          stepnow, ' (', trim(current_timestr) ,')', 'Next observation at time step', stepnow + nsteps, &
          ' (', trim(future_timestr) ,')'
          write(*,*) 'available | observation '
          do i=1,size(observations_assim)
            write(*,*) observations_assim(i)%ptr%is_available_at_next_assim_step, observations_assim(i)%ptr%name
          end do
    end if
  else
    ! *** End of assimilation process ***
     nsteps = 0
     doexit = 1

     IF (rank_world == 0) WRITE (*, '(i7, 3x, a)') &
          stepnow, 'No more observations - end assimilation'

  end if

  deallocate(obs_available)

END SUBROUTINE next_observation_pdaf
