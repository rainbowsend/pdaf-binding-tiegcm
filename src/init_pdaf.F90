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
! Top-level PDAF setup: configures and calls PDAF_init, wiring TIE-GCM state/config into the filter.

!BOP
!
! !ROUTINE: init_pdaf - Interface routine to call initialization of PDAF
!
! !INTERFACE:
SUBROUTINE init_pdaf()

! !DESCRIPTION:
! This routine collects the initialization of variables for PDAF.
! In addition, the initialization routine PDAF_init is called
! such that the internal initialization of PDAF is performed.
! This variant is for the online mode of PDAF.
!
! This routine is generic. However, it assumes a constant observation
! error (rms_obs). Further, with parallelization the local state
! dimension dimState is used.
!

  ! tie-gcm
  use hist_module, only:modeltime
  use input_module, only: pristart
  use mpi_module, only: TIEGCM_WORLD, mytid

  ! pdaf
  use pdaf, only: pdaf_init, pdaf_get_localfilter, pdaf_get_state, pdaf_set_debug_flag, &
                  PDAFomi_set_searchtype

  ! intern
  use configuration,&
      only: cfg_output, cfg_filter, cfg_ensemble, &
            assimilated_field_names, cfg_calibration
  use coordinates_module, only: init_coordinates_module
  use cell_id_coordinate_system, only: init_cell_id_coordinate_system
  use ensemble_module,&
      only: fill_state_0, init_ens_from_first_state
  use quantity_computation_module, only: init_quantity_computation_module
  use model_parameter_handling_module,&
    only: calibrated_parameter_names
  use model_parameter_IO_module,&
    only: model_parameter_writer
  use mod_assimilation,&
      only: screen, filtertype, subtype, dim_ens, &
             rms_obs, incremental, covartype, type_forget, forget, &
             rank_analysis_enkf, locweight, cutoff_radius, support_radius, &
             filename, type_trans, type_sqrt, delt_obs
  use mod_parallel_pdaf,&
      only: rank_world, n_modeltasks, task_id, &
            COMM_filter, COMM_couple, filterpe, abort_parallel
  use observations_module,&
      only: init_observations_module, observations
  use tiegcm_optimized_interpolator,&
      only: init_tiegcm_optimized_interpolator
  use quantity_info_module,&
      only: init_quantity_info_module
  use result_file_writer_module,&
      only: result_file_writer_module_set_precision, test_result_file_writer_module
  use result_writer_frontend,&
      only: init_result_writer_frontend, init_result_writer
  use state_module,&
      only: init_field_indices, state_vector,&
            idx_intern, nlevX, init_horizontal_coords_vector
  use structured_gird_subdomain_module,&
      only: init_structured_gird_subdomain_module, &
            disjoint_domain_regulary, &
            compute_mapping, &
            test_get_center_coords
  use time_module,&
      only: seconds_to_tgcm_steps, init_time_module

  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: tgcm.F
! Calls: init_pdaf_parse
! Calls: init_pdaf_info
! Calls: PDAF_init
! Calls: PDAF_get_state
!EOP

! Local variables
  INTEGER :: filter_param_i(7) ! Integer parameter array for filter
  REAL    :: filter_param_r(2) ! Real parameter array for filter
  INTEGER :: status_pdaf       ! PDAF status flag
  INTEGER :: doexit, steps     ! Not used in this implementation
  REAL    :: timenow           ! Not used in this implementation

  integer :: localfilter

  integer :: i

  ! External subroutines
  EXTERNAL :: init_ens         ! Ensemble initialization

  EXTERNAL :: next_observation_pdaf, & ! Provide time step, model time, 
                                       ! and dimension of next observation
       distribute_state_pdaf, &        ! Routine to distribute a state vector to model fields
       prepoststep_ens_pdaf            ! User supplied pre/poststep routine

! ***************************
! ***   Initialize PDAF   ***
! ***************************

  IF (rank_world == 0) THEN
     WRITE (*,'(/1x,a)') 'INITIALIZE PDAF - ONLINE MODE'
  END IF

  if( (n_modeltasks .lt. 2) .and. (cfg_filter%open_loop .eqv. .false.) ) then
    write(*,'(a,/,a,i3,a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!', &
                       'WARNING: ensemble size ( number of parallel model tasks ) is '&
                     //'smaller than two (', &
                     n_modeltasks, &
                     '). Prorgram is going to crash. Increase number or '&
                     //'set open_loop = .true.'
  end if

  CALL PDAF_set_debug_flag(0)

!   CALL PDAFomi_set_searchtype(0, 0)

  ! modeltime is initalized in advance.F Too late for next obs
  ! --> init here (same is done in advance again)
  modeltime(1:4) = pristart(:,1)

  call init_time_module()
  call init_coordinates_module()

  call init_field_indices()
  call init_horizontal_coords_vector()

  call init_quantity_info_module()

  if(cfg_calibration%apply) &
    call model_parameter_writer%create('results_calibration_parameters.nc')

  if(cfg_filter%open_loop .eqv. .false.) then
    ! *** State Options
    if(cfg_calibration%apply) then
      call state_vector%init( assimilated_field_names, &
                              name='state', &
                              calibrations=calibrated_parameter_names,&
                              print_map=.true. )
    else
      call state_vector%init( assimilated_field_names, &
                              name='state', &
                              print_map=.true. )
    end if
  end if

  call init_tiegcm_optimized_interpolator
  call init_cell_id_coordinate_system()
  call set_pdaf_omi_domain_limits()

  if(cfg_filter%open_loop .eqv. .false.) then

    ! *** Ensemble generation options
    init_ens_from_first_state = cfg_ensemble%overwrite_source
    if(cfg_ensemble%overwrite_source) then
      call fill_state_0()
    end if

! **********************************************************
! ***   CONTROL OF PDAF - used in call to PDAF_init      ***
! **********************************************************

  ! *** IO options ***
    screen      = 2  ! Write screen output (1) for output, (2) add timings

  ! *** Filter specific variables
  ! Type of filter
  !   (1) SEIK
  !   (2) EnKF
  !   (3) LSEIK
  !   (4) ETKF
  !   (5) LETKF
  !   (6) ESTKF
  !   (7) LESTKF
  !   (8) localized EnKF
  !   (9) NETF
  !  (10) LNETF
    filtertype = cfg_filter%filtertype

  ! Size of ensemble for all ensemble filters
  ! Number of EOFs to be used for SEEK
    dim_ens = n_modeltasks

  ! subtype of filter: 
  !   ESTKF:
  !     (0) Standard form of ESTKF
  !   LESTKF:
  !     (0) Standard form of LESTKF
    subtype = cfg_filter%subtype

  ! Type of ensemble transformation
  !   SEIK/LSEIK and ESTKF/LESTKF:
  !     (0) use deterministic omega
  !     (1) use random orthonormal omega orthogonal to (1,...,1)^T
  !     (2) use product of (0) with random orthonormal matrix with
  !         eigenvector (1,...,1)^T
  !   ETKF/LETKF:
  !     (0) use deterministic symmetric transformation
  !     (2) use product of (0) with random orthonormal matrix with
  !         eigenvector (1,...,1)^T
    type_trans = cfg_filter%type_trans

  ! Type of forgetting factor in SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
  !   (0) fixed
  !   (1) global adaptive
  !   (2) local adaptive for LSEIK/LETKF/LESTKF
    type_forget = cfg_filter%type_forget

  ! Forgetting factor
    forget  = cfg_filter%forget

  ! Type of transform matrix square-root
  !   (0) symmetric square root, (1) Cholesky decomposition
    type_sqrt = cfg_filter%type_sqrt

  ! (1) to perform incremental updating (only in SEIK/LSEIK!)
    incremental = 0

  ! Definition of factor in covar. matrix used in SEIK
  !   (0) for dim_ens^-1 (old SEIK)
  !   (1) for (dim_ens-1)^-1 (real ensemble covariance matrix)
  !   This parameter has also to be set internally in PDAF_init.
    covartype = 1

  ! rank to be considered for inversion of HPH
  ! in analysis of EnKF; (0) for analysis w/o eigendecomposition
    rank_analysis_enkf = cfg_filter%rank_analysis_enkf

! *********************************************************************
! ***   Settings for analysis steps  - used in call-back routines   ***
! *********************************************************************

  ! *** Forecast length (time interval between analysis steps) ***
  ! ATTENTION this variable has no impact and is only printed in init_pdaf_info
    delt_obs = seconds_to_tgcm_steps(cfg_filter%forecast_duration_sec)     ! Number of time steps between analysis/assimilation steps

  ! *** specifications for observations ***
    rms_obs = 0.5    ! Observation error standard deviation
                    ! for the Gaussian distribution 
  ! *** Localization settings
    locweight = cfg_filter%locweight ! Type of localizating weighting
                      ! see subroutine PDAFomi_weights_l in PDAFomi_obs_l.F90
                      ! and subroutine PDAF_local_weight in PDAF_local_weight.F90
                      !   (0) constant weight of 1
                      !   (1) exponentially decreasing with support_radius
                      !       e.g. support_radius = H/ln(2) -> H: half-life
                      !   (2) use 5th-order polynomial
                      !   (3) regulated localization of R with mean error variance
                      !   (4) regulated localization of R with single-point error variance
    cutoff_radius = cfg_filter%cutoff_radius
    !
    ! delta model height ~ 400 km
    ! @ 5.0 deg 100 km altiude at equator -> cell size = (6371 km+100 km)*5/180*pi = 564 km
    ! @ 2.5 deg 100 km altiude at equator -> cell size = (6371 km+100 km)*5/180*pi = 282 km
    !
    ! cutoff_radius = 6371E+3*30/180*3.14  ! Range in meters
    support_radius = cfg_filter%support_radius  ! Support range for 5th-order polynomial
                            ! or range for 1/e for exponential weighting

  ! *** File names
    filename = 'output.dat'


! ***********************************
! *** Some optional functionality ***
! ***********************************

  ! *** Parse command line options   ***
  ! *** This is optional, but useful ***

  ! this overwrites support_radius and is not needed since we use namelist files
  !  call init_pdaf_parse()

  ! *** Initial Screen output ***
  ! *** This is optional      ***

    IF (rank_world == 0) call init_pdaf_info()


! *****************************************************
! *** Call PDAF initialization routine on all PEs.  ***
! ***                                               ***
! *** Here, the full selection of filters is        ***
! *** implemented. In a real implementation, one    ***
! *** reduce this to selected filters.              ***
! ***                                               ***
! *** For all filters, first the arrays of integer  ***
! *** and real number parameters are initialized.   ***
! *** Subsequently, PDAF_init is called.            ***
! *****************************************************

    whichinit: IF (filtertype == 2) THEN
      ! *** EnKF with Monte Carlo init ***
      filter_param_i(1) = state_vector%map%size(rank=mytid) ! State dimension
      filter_param_i(2) = dim_ens       ! Size of ensemble
      filter_param_i(3) = rank_analysis_enkf ! Rank of speudo-inverse in analysis
      filter_param_i(4) = incremental ! Whether to perform incremental analysis
      filter_param_i(5) = 0           ! Smoother lag (not implemented here)
      filter_param_r(1) = forget      ! Forgetting factor

      CALL PDAF_init(filtertype, subtype, 0, &
            filter_param_i, 6,&
            filter_param_r, 2, &
            TIEGCM_WORLD%MPI_VAL, COMM_filter%MPI_VAL, COMM_couple%MPI_VAL, &
            task_id, n_modeltasks, filterpe, init_ens, &
            screen, status_pdaf)
    ELSE
      ! *** All other filters                       ***
      ! *** SEIK, LSEIK, ETKF, LETKF, ESTKF, LESTKF ***
      filter_param_i(1) = state_vector%map%size(rank=mytid) ! State dimension
      filter_param_i(2) = dim_ens     ! Size of ensemble
      filter_param_i(3) = 0           ! Smoother lag (not implemented here)
      filter_param_i(4) = incremental ! Whether to perform incremental analysis (SEIK only)
      filter_param_i(5) = type_forget ! Type of forgetting factor
      filter_param_i(6) = type_trans  ! Type of ensemble transformation
      filter_param_i(7) = type_sqrt   ! Type of transform square-root (SEIK-sub4/ESTKF)
      filter_param_r(1) = forget      ! Forgetting factor

      CALL PDAF_init(filtertype, subtype, 0, &
            filter_param_i, 7,&
            filter_param_r, 2, &
            TIEGCM_WORLD%MPI_VAL, COMM_filter%MPI_VAL, COMM_couple%MPI_VAL, &
            task_id, n_modeltasks, filterpe, init_ens, &
            screen, status_pdaf)
    END IF whichinit

  ! *** Check whether initialization of PDAF was successful ***
    IF (status_pdaf /= 0) THEN
      WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
            'ERROR ', status_pdaf, &
            ' in initialization of PDAF - stopping! (PE ', rank_world,')'
      CALL abort_parallel()
    END IF

  end if

  ! ******************************'***
  ! Domain localization of model grid
  ! ******************************'***
  CALL PDAF_get_localfilter(localfilter)
  if((cfg_filter%open_loop .eqv. .false.) .and. (localfilter==1)) then
      ! localization
    call init_structured_gird_subdomain_module(nlevX,&
                                               idx_intern(mytid)%nlons,&
                                               idx_intern(mytid)%nlats)
    ! lev lon lat
    call disjoint_domain_regulary(cfg_filter%sub_domain_size_vertical,&
                                  cfg_filter%sub_domain_size_zonal,&
                                  cfg_filter%sub_domain_size_meridional)
    call compute_mapping
!     call test_get_center_coords
  end if

! ******************************'***
! *** Prepare observations ***
! ******************************'*** 
! must be called before PDAF_get_state, which calls next_observation_pdaf
! next_observation_pdaf requires info read by init_observation_settings 

  call result_file_writer_module_set_precision(cfg_output%use_double_precision)
  call init_result_writer_frontend()
!   call test_result_file_writer_module

  call init_observations_module()

  if(cfg_filter%open_loop .eqv. .false.) then

    ! observation type table
    write(*,*) 'pdaf omi observation type table'
    write(*,*) 'type        used'
    if(allocated(observations))then
      do i=1,size(observations)
        write(*,'(a12,l2)') trim(observations(i)%ptr%name), observations(i)%ptr%assimilate
      end do
    end if
  end if

  ! must be called after initalization routines for omi
  call init_quantity_computation_module

  ! must be called after init_result_writer_frontend and initalization routines for omi
  call init_result_writer
  call allocate_diagnostic_fields

! ******************************'***
! *** Prepare ensemble forecasts ***
! ******************************'***

if( cfg_filter%open_loop .eqv. .false. ) then
 CALL PDAF_get_state(steps, timenow, doexit, next_observation_pdaf, &
      distribute_state_pdaf, prepoststep_ens_pdaf, status_pdaf)
end if

END SUBROUTINE init_pdaf

! due to circle dependency this routine is located here
subroutine allocate_diagnostic_fields

  ! tie-gcm
  use fields_module, only: levd0,levd1,lond0,lond1,latd0,latd1

  ! intern
  use quantity_computation_module, only: qset
  use quantity_info_module, only: additional_fields, additional_fields_2d, quantity_infos, LEVEL_NONE

  implicit none

  character(len=16), dimension(:), allocatable :: diagnostic_fields
  integer :: n
  integer :: n_3d, n_2d
  integer :: i
  integer :: i_3d, i_2d
  integer :: fidx


  call qset%diagnostic%to_array(diagnostic_fields)
  n = size(diagnostic_fields,dim=1)

  ! quantities without vertical extent are stored with a single level only
  n_2d = 0
  do i =1,n
    fidx = findloc(quantity_infos%name,diagnostic_fields(i),dim=1)
    if(quantity_infos(fidx)%level==LEVEL_NONE) n_2d = n_2d + 1
  end do
  n_3d = n - n_2d

  allocate(additional_fields(levd0:levd1,lond0:lond1,latd0:latd1,n_3d))
  allocate(additional_fields_2d(levd0:levd0,lond0:lond1,latd0:latd1,n_2d))

  i_3d = 0
  i_2d = 0
  do i =1,n
    fidx = findloc(quantity_infos%name,diagnostic_fields(i),dim=1)
    ! ATTENTION explicitly initalize bounds of the pointer. Some compiler version default them to 1:size
    if(quantity_infos(fidx)%level==LEVEL_NONE) then
      i_2d = i_2d + 1
      quantity_infos(fidx)%data(levd0:levd0,lond0:lond1,latd0:latd1) => additional_fields_2d(:,:,:,i_2d)
    else
      i_3d = i_3d + 1
      quantity_infos(fidx)%data(levd0:levd1,lond0:lond1,latd0:latd1) => additional_fields(:,:,:,i_3d)
    end if
  end do

  if(allocated(diagnostic_fields)) deallocate(diagnostic_fields)

end subroutine
