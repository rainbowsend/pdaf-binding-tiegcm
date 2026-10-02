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
! Builds the PE-local PDAF state vector from TIE-GCM fields and precomputes interpolated observation equivalents.

!BOP
!
! !ROUTINE: collect_state_pdaf --- Initialize state vector from model fields
!
! !INTERFACE:
SUBROUTINE collect_state_pdaf(dim_p, state_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! This subroutine is called during the forecast 
! phase from PDAF\_put\_state\_X or PDAF\_assimilate\_X
! after the propagation of each ensemble member. 
! The supplied state vector has to be initialized
! from the model fields (typically via a module). 
! With parallelization, MPI communication might be 
! required to initialize state vectors for all 
! subdomains on the model PEs. 
!
! The routine is executed by each process that is
! participating in the model integrations.
!
! ATTENTION this routine is not only collecting the current state, but also
! interpolating the state to the observations. Actually, pdaf intends the 
! 'ops_op_pdaf' routine for this. However, ops_op_pdaf is only called by filter
! processes. Since the interpolation weights are different at each ensemble
! member obs_op would have to calculate the weights of all members 
! sequentially. To speed things up we calculate everything in this routine--
! which is called by all processes--and gather the results on the filter 
! processes
!
! !USES:

  ! extern
  use mpi_f08

  ! tie-gcm
  use fields_module, only: itc
  use mpi_module, only: mytid

  ! intern
  use array_print_module,&
      only: printMat
  use configuration,&
      only: cfg_calibration
  use georeferenced_data_module,  only: map_bundle_to_vec
  use mod_assimilation, only:&
      analysis_step_count
  use mod_parallel_pdaf,&
      only: COMM_couple
  use observations_module,&
      only: msis_cal, tum_ne, tum_vtec, satellite
  use quantity_computation_module,&
      only: compute_results
  use result_writer_frontend,&
      only: compute_and_write_results
  use state_module,&
      only: state_vector
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p           ! PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  ! local state vector

  ! local
  integer :: ierr

  real, pointer :: den

  real, allocatable, dimension(:) :: m_state_p

  integer :: dim_obs_p

  integer :: i

! !CALLING SEQUENCE:
! Called by: PDAF_put_state_X    (as U_coll_state)
! Called by: PDAF_assimilate_X   (as U_coll_state)
!EOP

! *************************************************
! *** Initialize state vector from model fields ***
! *************************************************

  write(*,*) "Initialize state vector from model (", analysis_step_count,")"

  call state_vector%fill_state_p(state_p, itc)

  if(cfg_calibration%apply)then
    call state_vector%fill_state_p_calibration(state_p)
  end if

  ! Compute observations (and other quantities that are written to result file)
  call compute_results(.true.)

  ! First time (analysis_step_count==0) is not for performaing analysis
  if(analysis_step_count>0)then
    call compute_and_write_results("forecast")
  end if


! ******************************************************
! *** perform interpolation of state to observations ***
! ******************************************************
! PDAF calls the observation operator sequentially for all members on the filter processes.
! To speed up the computation we calculate the interpolation (of msis) here and send the results to
! the filter processes.

! ATTENTION
! The result of the observation operator are "computed observations"
! The observation operator consists of two steps
! (1) computation of observed quantity from state (if state is not observed directly)
! (2) interpolation (if observations and state are not located at same time and position)
! Since interpolation is not exact, changing the order of the steps does affect the results
! We do the interpolation always at the end, so that the computation of quantities is not 
! affected by interpolation errors

! The observations are computed in
! call compute_results(.true.)
! Here we only need to collect the result

  do i=1,size(msis_cal)
    if(msis_cal(i)%assimilate)then
      write(*,*) "observation operator for ", msis_cal(i)%config%name

      dim_obs_p = msis_cal(i)%tme_obs%map%size_R(mytid)
      allocate( m_state_p( dim_obs_p ) )

      call map_bundle_to_vec(msis_cal(i)%tme_obs, msis_cal(i)%dst%data, m_state_p)

      ! Gather interpolated values at filter
      ! computed_msis_observations: local observation vector size X ensemble size
      call MPI_Gather( m_state_p,  dim_obs_p, MPI_DOUBLE, &
                      msis_cal(i)%computed_observations, dim_obs_p, MPI_DOUBLE, &
                      0, &
                      COMM_couple, ierr)
      deallocate( m_state_p )

    end if
  end do

  if(tum_ne%assimilate)then
    write(*,*) "observation operator for ", tum_ne%config%name

    dim_obs_p = tum_ne%tme_obs%map%size_R(mytid)
    allocate( m_state_p( dim_obs_p ) )

    call map_bundle_to_vec(tum_ne%tme_obs, tum_ne%dst%data, m_state_p)

    ! Gather interpolated values at filter
    ! computed_msis_observations: local observation vector size X ensemble size
    call MPI_Gather( m_state_p,  dim_obs_p, MPI_DOUBLE, &
                    tum_ne%computed_observations, dim_obs_p, MPI_DOUBLE, &
                    0, &
                    COMM_couple, ierr)
    deallocate( m_state_p )

  end if

  if(tum_vtec%assimilate)then
    write(*,*) "observation operator for ", tum_vtec%config%name

    dim_obs_p = tum_vtec%tme_obs%map%size_R(mytid)
    allocate( m_state_p( dim_obs_p ) )

    call map_bundle_to_vec(tum_vtec%tme_obs, tum_vtec%dst%data, m_state_p)

    ! Gather interpolated values at filter
    ! computed_msis_observations: local observation vector size X ensemble size
    call MPI_Gather( m_state_p,  dim_obs_p, MPI_DOUBLE, &
                    tum_vtec%computed_observations, dim_obs_p, MPI_DOUBLE, &
                    0, &
                    COMM_couple, ierr)
    deallocate( m_state_p )

  end if

  do i=1,size(satellite)
    if(satellite(i)%assimilate)then
      write(*,*) "observation operator for ", satellite(i)%config%name

      dim_obs_p = 1
      allocate( m_state_p( dim_obs_p ) )

      den => null()
      ! den is available at all ranks (WORLD)
      call satellite(i)%dst%data%get(satellite(i)%observed_quantity,den)

      m_state_p(1) = den

      ! Gather interpolated values at filter
      call MPI_Gather( m_state_p,  dim_obs_p, MPI_DOUBLE, &
                      satellite(i)%computed_observations, dim_obs_p, MPI_DOUBLE, &
                      0, &
                      COMM_couple, ierr)
      deallocate( m_state_p )

    end if
  end do

END SUBROUTINE collect_state_pdaf
