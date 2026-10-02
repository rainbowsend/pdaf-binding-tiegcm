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
! PDAF-OMI observation type for along-track satellite density observations, including trajectory interpolation and localization.

!
!
!                      -------------------------------
!                      |  observation_interface      |
!                      |-----------------------------|
!                      |  TYPE(obs_f) :: full_obs    |
!                      |  integer :: obs_l_id        |
!                      |  LOGICAL :: assimilate      |
!                      |  REAL    :: rms_obs         |
!                      |  character(len=64) :: name  |
!                      |-----------------------------|
!                      |  init_dim_obs               |
!                      |  obs_op                     |
!                      |  init_dim_obs_l             |
!                      |  localize_covar             |
!                      -------------------------------
!                                    |
!                                    |
!                                    V
!     --------------------------------------------------------------
!     |                 tgcm_pdaf_omi_obs_type_module              |
!     |------------------------------------------------------------|
!     |real, dimension(:,:), allocatable :: computed_observations  |
!     |----------------------------------------------------------- |
!     | next_observation                                           |
!     | init_full_observation_dimension                            |
!     --------------------------------------------------------------
!                          |                      |
!                          |                      |
!                          V                      V
!    +++++++++++++++++++++++++++++++         -------------------------------
!    +       obs_satellite         +         |          obs_tme_grid       |
!    +-----------------------------+         |-----------------------------|
!    + type(trajectory_data)::data +         |        ...                  |
!    +-----------------------------+         |-----------------------------|
!    +                             +         |        ...                  |
!    +        ...                  +         -------------------------------
!    +++++++++++++++++++++++++++++++            |                      |   |___________________________________
!                                               |                      |                                       |
!                                               V                      V                                       V
!                         -------------------------------      -------------------------------       -------------------------------
!                         |          obs_den_grid       |      |          obs_tum_ne         |       |          obs_tum_vtec       |
!                         |-----------------------------|      |-----------------------------|       |-----------------------------|
!                         |        ...                  |      |        ...                  |       |        ...                  |
!                         -------------------------------      -------------------------------       -------------------------------
!
MODULE obs_satellite_pdafomi

  ! intern
  use trajectory_data_module, only: trajectory_data
  use configuration, only: satellite_file
  use georeferenced_data_module, only: point
  use result_file_writer_module, only: nc_trajectory
  use tgcm_pdaf_omi_obs_type_module, only: tgcm_pdaf_observation_interface

  ! tie-gcm
  use fields_module, only: shortname_len

  implicit none

  type, extends(tgcm_pdaf_observation_interface) :: obs_satellite
    type(trajectory_data) :: data
    type(satellite_file) :: config
    real :: time
    real, dimension(3) :: position ! lon(-180 deg : 180 deg), lat(-90 deg : 90 deg), alt
    logical :: obs_on_this_rank
    type(nc_trajectory) :: spatial_domain
    type(point) :: dst
    character(len=shortname_len) :: observed_quantity
    contains
    ! derived procedures
    procedure, pass(this) :: init_dim_obs_l => obs_satellite_init_dim_obs_l
    procedure, pass(this) :: localize_covar => obs_satellite_localize_covar
    procedure, pass(this) :: next_observation => obs_satellite_next_observation
    procedure, pass(this) :: obs_op => obs_satellite_obs_op
    procedure, pass(this) :: init_dim_obs => obs_satellite_init_dim_obs
    procedure, pass(this) :: deallocate => obs_satellite_deallocate
    ! class procedures
    procedure, pass(this) :: init => obs_satellite_init
  end type

  contains

  !> Loads the satellite trajectory/density file (format from config), sets up its result-writer domain, and allocates per-member computed-observation storage if this type is assimilated or always saved.
  subroutine obs_satellite_init(this, config)

    ! intern
    use configuration, only: satellite_file, cfg_output
    use mod_parallel_pdaf, only: n_modeltasks
    use tgcm_pdaf_omi_obs_type_module, only: get_tiegcm_native_grid

    implicit none

    ! arguments
    class(obs_satellite) :: this
    type(satellite_file), intent(in) :: config

    this%config=config
    this%name = config%name

    this%assimilate = config%apply

    this%full_obs%use_global_obs = 0

    select case(config%file_format)
      case('igg','IGG')
          call this%data%read_denswind(config%file)
          this%observed_quantity = 'DEN'
      case('toleos','toleos_reduced','toleos_short','groops')
          call this%data%read_txt(config%file,config%file_format)
          this%observed_quantity = 'DEN'
      case default
        call shutdown(config%file_format//' is an invalid file format for satellite observation')
    end select

    this%tiegcm_native_grid = get_tiegcm_native_grid(this%observed_quantity)


    call this%spatial_domain%init(domain_name=trim(config%name),&
                                  save_members=cfg_output%save_members,&
                                  kmax=cfg_output%max_moment,&
                                  write_every_sec=config%write_every_sec,&
                                  force_write_on_update=config%force_write_on_update,&
                                  save_n_steps_after_update=cfg_output%save_n_steps_after_update)
    call this%spatial_domain%link_dataset(this%data)

    ! TODO rename dst
    call this%dst%init(this%data,this%spatial_domain,this%name)

    call this%link_to_writer(this%spatial_domain,  (/this%observed_quantity/))

    if ((this%assimilate) .or. (config%always_save)) then
      ! the interpolation is already done in collect state, so that all processes can do it instead of
      ! sequential execution here on filter prcocesses
      allocate( this%computed_observations(1, n_modeltasks ) )
    end if

  end subroutine

  !> Deallocates this observation's trajectory data, computed-observations buffer, and output-mapping resources.
  subroutine obs_satellite_deallocate(this)

    implicit none

    ! arguments
    class(obs_satellite) :: this

    call this%data%deallocate
    if (allocated (this%computed_observations))   deallocate (this%computed_observations)

    call this%dst%deallocate()

  end subroutine

  !> Reports whether this satellite type is assimilated and has a trajectory sample valid at the given next analysis time.
  function obs_satellite_next_observation(this, next_analysis_step) result(is_available)

    ! extern
    use esmf

    ! intern
    use time_module, only: convert_to_seconds_since_ref

    implicit none
    ! arguments
    class(obs_satellite) :: this
    type(ESMF_Time), intent(in) :: next_analysis_step

    ! result
    logical :: is_available

    is_available = .false.
    if( this%assimilate .eqv. .true. ) then
      call convert_to_seconds_since_ref(next_analysis_step, this%time)
      is_available = this%data%valid_at_epoch(this%time)
    end if
  end function

  !> Sets the local observation dimension for specified domain.
  SUBROUTINE obs_satellite_init_dim_obs_l(this, domain_p, step, dim_obs, dim_obs_l)

    ! PDAF
    use PDAF, only: PDAFomi_init_dim_obs_l
    use PDAFomi_obs_l, only: obs_l

    ! intern
    use mod_assimilation, only: coords_l, cutoff_radius, locweight, support_radius

    IMPLICIT NONE

! *** Arguments ***
    class(obs_satellite) :: this
    INTEGER, INTENT(in)  :: domain_p     !< Index of current local analysis domain
    INTEGER, INTENT(in)  :: step         !< Current time step
    INTEGER, INTENT(in)  :: dim_obs      !< Full dimension of observation vector
    INTEGER, INTENT(inout) :: dim_obs_l  !< Local dimension of observation vector

! *** local variables ***
    type(obs_l), pointer :: local_obs    ! this thread's local observation


!     ! Template reminder - delete when implementing functionality
!     WRITE (*,*) 'TEMPLATE init_msis_pdafomi_TEMPLATE.F90: Initialize local observations'

! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

    ! Here one has to specify the coordinates of the local analysis domain
    ! (coords_l) and the localization variables, which can be different for
    ! each observation type and can be made dependent on the index DOMAIN_P.
    ! coords_l should be set in the call-back routine init_dim_l.

    call this%local_obs(local_obs)

    ! coords_l and the radii are always given for the full 3d case. They are
    ! restricted to ncoord, which is 2 for observations without vertical
    ! extent (case `LEVEL_NONE`)
    CALL PDAFomi_init_dim_obs_l(local_obs, this%full_obs, &
         coords_l(1:this%full_obs%ncoord), locweight, &
         cutoff_radius(1:this%full_obs%ncoord), &
         support_radius(1:this%full_obs%ncoord), dim_obs_l)

  END SUBROUTINE obs_satellite_init_dim_obs_l

  !> Applies covariance localization.
  SUBROUTINE obs_satellite_localize_covar(this, dim_p, dim_obs, HP_p, HPH, coords_p)

    ! PDAF
    use PDAF, only: PDAFomi_localize_covar

    ! intern
    use mod_assimilation, only: cutoff_radius, locweight, support_radius

    IMPLICIT NONE

! *** Arguments ***
    class(obs_satellite) :: this
    INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
    INTEGER, INTENT(in) :: dim_obs               !< Dimension of observation vector
    REAL, INTENT(inout) :: HP_p(dim_obs, dim_p)  !< PE local part of matrix HP
    REAL, INTENT(inout) :: HPH(dim_obs, dim_obs) !< Matrix HPH
    REAL, INTENT(in)    :: coords_p(:,:)         !< Coordinates of state vector elements


    ! Template reminder - delete when implementing functionality
    WRITE (*,*) 'TEMPLATE init_msis_pdafomi_TEMPLATE.F90: Apply covariance localization'

! *************************************
! *** Apply covariance localization ***
! *************************************

    ! Here one has to specify the three localization variables
    ! which can be different for each observation type.

    CALL PDAFomi_localize_covar(this%full_obs, dim_p, locweight, &
         cutoff_radius(1:this%full_obs%ncoord), &
         support_radius(1:this%full_obs%ncoord), &
         coords_p, HP_p, HPH)

  END SUBROUTINE obs_satellite_localize_covar

  !> Applies the observation operator: computes this ensemble member's satellite density observation equivalent (using precomputed values for non-mean members) and gathers it into the full observed-state vector.
  SUBROUTINE obs_satellite_obs_op(this, dim_p, dim_obs, state_p, ostate)

    ! pdaf
    use PDAFomi_obs_f, only:  PDAFomi_gather_obsstate, debug
    use PDAF, only: PDAF_get_obsmemberid

    ! intern
    use quantity_computation_module, only: zg_mid_mean
    use state_module, only:  state_vector


    IMPLICIT NONE

! *** Arguments ***
    class(obs_satellite) :: this
    INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
    INTEGER, INTENT(in) :: dim_obs               !< Dimension of full observed state (all observed fields)
    REAL, INTENT(in)    :: state_p(dim_p)        !< PE-local model state
    REAL, INTENT(inout) :: ostate(dim_obs)       !< Full observed state

! local

    REAL, ALLOCATABLE :: ostate_p(:)       ! local observed part of state vector
    integer :: current_model_instance

    IF (debug>0) THEN
      WRITE (*,*) '++ OMI-debug: ', debug, 'obs_op', this%name,'-- START'
    END IF

! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

    IF (this%full_obs%doassim==1) THEN

        call PDAF_get_obsmemberid( current_model_instance )

        write(*,*) "observation operator for ensemble member ", current_model_instance

        ALLOCATE(ostate_p(this%full_obs%dim_obs_p))

        if(current_model_instance == 0) then
          ! apply obs op to mean state (only global filters)

          ! This must be called from all ranks for halo exchange, however observation is only on one rank. Thus, this is only called at the rank including the obs with original PDAF
          call obs_satellite_obs_op_computation(this, state_p, state_vector, ostate_p(1), zg_mid=zg_mid_mean)

        else
          ! computed in collect_state_pdaf
           ostate_p  = this%computed_observations(1, current_model_instance)
        end if

       ! *************************************************
       ! *** Global: Gather full observed state vector ***
       ! ***            THIS IS MANDATORY!             ***
       ! *************************************************
       CALL PDAFomi_gather_obsstate(this%full_obs, ostate_p, ostate)

        DEALLOCATE(ostate_p)
    END IF

    ! Print debug information
    IF (debug>0) &
      WRITE (*,*) '++ OMI-debug: ', debug, 'obs_op_msis -- END'

  END SUBROUTINE obs_satellite_obs_op

  !> Computes the density (DEN) observation equivalent at this satellite's position by interpolating the given state onto its location.
  subroutine obs_satellite_obs_op_computation(this, state_p, state_map, computed_obs, zg_mid, zg_mid_nm, zg_int, zg_int_nm)

    ! tie-gcm
    use fields_module,only: levd0,levd1,lond0,lond1,latd0,latd1

    ! intern
    use configuration, only: cfg_filter
    use quantity_info_module, only: quantity_info, get_info
    use state_module, only:  state_vector_mapping
    use tiegcm_optimized_interpolator, only: sparse_state_interpolator

    implicit none

    ! arguments
    class(obs_satellite) :: this
    real, dimension(:), intent(in) :: state_p
    type(state_vector_mapping), intent(in) :: state_map
    real, intent(out) :: computed_obs
    real, dimension(levd0:levd1,lond0:lond1,latd0:latd1), intent(inout), optional :: zg_mid, zg_mid_nm, zg_int, zg_int_nm

    ! local
    type(quantity_info), pointer :: quantity
    type(sparse_state_interpolator) :: interp

    quantity=>get_info(this%observed_quantity)
    call quantity%calc(state_map, state_p)

    call interp%init(this%position,  zg_mid=zg_mid, degree=cfg_filter%spline_degree )
    call interp%interpolate(quantity,quantity%data,computed_obs)

    if(this%obs_on_this_rank) then
      write(*,*) 'computed_obs: ', computed_obs
    end if

    call interp%deallocate()

  end subroutine obs_satellite_obs_op_computation


  !> Reads this satellite's observation value at the current time, determines its owning subdomain and coordinates, sets its observation error, and gathers the full PDAF-OMI observation arrays.
  SUBROUTINE obs_satellite_init_dim_obs(this, step, dim_obs)

    ! extern
    use esmf
    use mpi_f08

    ! tie-gcm
    use cons_module, only: pi

    ! pdaf
    use PDAF, only: PDAFomi_gather_obs

    ! intern
    use array_print_module, only: printMat
    use array_mapping_module, only: flatten, map_to_3d
    use cell_id_coordinate_system, only: to_cell_id_coordinate, to_zonal_meridional_vertical
    use configuration, only: cfg_filter
    use mod_assimilation, only: cutoff_radius
!     use mod_parallel_pdaf, only: rank_filter
    use tgcm_pdaf_omi_obs_type_module, only: COORD_SPH_3D, COORD_CELL_IDX_3D
    use tiegcm_optimized_interpolator, only: is_within_sub_domain, check_unambigous_sub_domain_assignment
    use quantity_computation_module, only: zg_mid_mean
    use quantity_info_module, only: LEVEL_MID

    IMPLICIT NONE

! *** Arguments ***
    class(obs_satellite) :: this
    INTEGER, INTENT(in)    :: step       !< Current time step
    INTEGER, INTENT(inout) :: dim_obs    !< Dimension of full observation vector

! *** Local variables ***
    INTEGER :: dim_obs_p                 ! Number of process-local observations
    REAL, ALLOCATABLE :: obs_p(:)        ! PE-local observation vector
    REAL, ALLOCATABLE :: ivar_obs_p(:)   ! PE-local inverse observation error variance
    REAL, ALLOCATABLE :: ocoord_p(:,:)   ! PE-local observation coordinates

    real, dimension(3) :: position_cell_id_coordinate
    real :: val

! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

    call this%init_full_observation_dimension()

! **********************************
! *** Read PE-local observations ***
! **********************************
    call this%data%get_at_epoch(this%time,val,this%position)

    ! convert kg/m3 to g/cm3
    val = val/1000

    write(*,*) 'time: ', this%time, ' pos: ', this%position , ' val: ', val

    this%obs_on_this_rank = is_within_sub_domain(lon=this%position(1),lat=this%position(2))

    ! call check_unambigous_sub_domain_assignment(this%obs_on_this_rank)

! ***********************************************************
! *** Count available observations for the process domain ***
! *** and initialize index and coordinate arrays.         ***
! ***********************************************************

    ! *** Count valid observations that lie within the process sub-domain ***

    if(this%obs_on_this_rank) then
      write(*,*) 'obs is in this sub domain'

      dim_obs_p = 1

      ! *** Initialize vector of observations on the process sub-domain ***

      ALLOCATE(obs_p(dim_obs_p))
      obs_p = val
    else
      dim_obs_p = 0
      allocate(obs_p(0))
    end if

    ! *** Initialize coordinate array of observations on the process sub-domain ***
    ! ATTENTION assume dim_obs_p==1
    ALLOCATE(ocoord_p(this%full_obs%ncoord, 1))
    select case(cfg_filter%localization_coord_sys)
      case(COORD_SPH_3D)
        ! lon lat alt
        ocoord_p(1,1) = this%position(1)/180.*pi ! radians east
        ocoord_p(2,1) = this%position(2)/180.*pi ! geocentric, radians north
        ocoord_p(3,1) = this%position(3)         ! meters above ellipsoid
      case(COORD_CELL_IDX_3D)
        ! ATTENTION calls interpolator that requires all model ranks!
        position_cell_id_coordinate = to_cell_id_coordinate(zg_mid_mean, LEVEL_MID, this%position)
        call to_zonal_meridional_vertical(position_cell_id_coordinate)
        ocoord_p(:,1) = position_cell_id_coordinate
    end select

!     write(*,*) 'ocoord_p: ', ocoord_p(:,1)

    ! not needed in that implementation
    ALLOCATE(this%full_obs%id_obs_p(0,0))

    if(this%obs_on_this_rank) then
  ! ****************************************************************
  ! *** Define observation errors for process-local observations ***
  ! ****************************************************************

      ALLOCATE(ivar_obs_p(dim_obs_p))

     ! std is 10% of value
     ivar_obs_p = val/10

     ! apply weight
     ivar_obs_p = ivar_obs_p/this%config%weight

     ! std to inverse variance
     ivar_obs_p = 1.0/(ivar_obs_p**2)

    else
      allocate(ivar_obs_p(0))
    end if

! ****************************************
! *** Gather global observation arrays ***
! ****************************************

    ! This routine is generic for the case that only the observations,
    ! inverse variances and observation coordinates are gathered

    CALL PDAFomi_gather_obs(this%full_obs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
         this%full_obs%ncoord, maxval(cutoff_radius(1:this%full_obs%ncoord)), dim_obs)

    ! *********************************************************
    ! *** For twin experiment: Read synthetic observations  ***
    ! *********************************************************

    !   IF (twin_experiment .AND. filtertype/=11) THEN
    !      CALL read_syn_obs(file_syntobs_TYPE, dim_obs, full_obs%obs_f, 0, 1-rank_filter)
    !   END IF


    ! ********************
    ! *** Finishing up ***
    ! ********************

    ! Deallocate all local arrays
    DEALLOCATE(obs_p, ocoord_p, ivar_obs_p)

    ! Arrays in THISOBS have to be deallocated after the analysis step
    ! by a call to deallocate_obs() in prepoststep_pdaf.

  END SUBROUTINE obs_satellite_init_dim_obs

end module obs_satellite_pdafomi
