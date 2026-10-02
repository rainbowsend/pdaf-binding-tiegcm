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
! PDAF-OMI observation types for gridded data: base regular-grid type plus MSIS-calibrated density and TUM electron-density observations.

!
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
!                                    |
!                                    |
!                                    V
!                      -------------------------------
!                      |          obs_tme_grid       |
!                      |-----------------------------|
!                      |        ...                  |
!                      -------------------------------
!                          |                      |  |_________________________________________
!                          |                      |                                            |
!                          V                      V                                            v
!    -------------------------------      -------------------------------     -------------------------------
!    |          obs_den_grid       |      |          obs_tum_ne         |     |          obs_tum_vtec       |
!    |-----------------------------|      |-----------------------------|     |-----------------------------|
!    |        ...                  |      |        ...                  |     |        ...                  |
!    -------------------------------      -------------------------------     -------------------------------
!
MODULE obs_tme_grid_pdafomi

  ! intern
  use configuration, only: key_len
  use georeferenced_data_module, only: regular_grid
  use grid_observation_module, only: reg_grid_dataset_root, reg_grid_dataset_group
  use tgcm_pdaf_omi_obs_type_module, only: tgcm_pdaf_observation_interface
  use result_file_writer_module, only: nc_reg_grid

  implicit none

  type, abstract, extends(tgcm_pdaf_observation_interface) :: obs_tme_grid
    type(reg_grid_dataset_root) :: tme_data_file
    type(reg_grid_dataset_group) :: tme_obs, tme_std
    type(nc_reg_grid) :: spatial_domain
    type(regular_grid) :: dst
    character(len=key_len), dimension(:), allocatable :: field_names
    contains
    ! derived procedures
    procedure, pass(this) :: init_dim_obs_l => obs_tme_grid_init_dim_obs_l
    procedure, pass(this) :: localize_covar => obs_tme_grid_localize_covar
    procedure, pass(this) :: next_observation => obs_tme_grid_next_observation
    procedure, pass(this) :: obs_op => obs_tme_grid_obs_op
    ! Implementation specific procedures
    procedure, pass(this) :: obs_op_computation => obs_tme_grid_obs_op_computation
    procedure, pass(this) :: init_coordinate_arrays => obs_tme_grid_init_coordinate_arrays
  end type

  contains

  !> Reports whether this grid is assimilated and has a sample valid at the given next analysis time.
  function obs_tme_grid_next_observation(this, next_analysis_step) result (is_available)

    ! extern
    use esmf

    ! intern
    use time_module, only: binary_search_esmf_time

    implicit none
    ! arguments
    class(obs_tme_grid) :: this
    type(ESMF_Time), intent(in) :: next_analysis_step

    ! result
    logical :: is_available

    ! local
    integer :: t_idx_prev
    integer :: n_epochs
    integer :: t_idx
    character(len=80) :: future_timestr

    is_available = .false.

    if( this%assimilate .eqv. .true. ) then

      call ESMF_TimeGet(next_analysis_step, timeString=future_timestr)

      t_idx_prev = this%tme_data_file%t_idx
      n_epochs = size(this%tme_data_file%epochs)
      t_idx = binary_search_esmf_time(this%tme_data_file%epochs(t_idx_prev:n_epochs),next_analysis_step)

      if(t_idx .gt. 0) then
        t_idx = t_idx_prev + t_idx -1
        this%tme_data_file%t_idx = t_idx
        is_available = .true.
        write(*,'(5a,i6)') 'next observation ', trim(this%name), ' ', trim(future_timestr), ' t_idx=', t_idx
      else
        write(*,*) 'no observation available for ',&
                              trim(this%name), &
                              '. Period is : ',&
                              trim(this%tme_data_file%first_timestr),&
                              ' - ',&
                              trim(this%tme_data_file%last_timestr), &
                              '   requested epoch is: ',&
                              trim(future_timestr)
        is_available = .false.
      end if
    end if
  end function

 !> Applies the observation operator to a state vector (using precomputed interpolated values for ensemble members) and gathers the full observed state.
 SUBROUTINE obs_tme_grid_obs_op(this, dim_p, dim_obs, state_p, ostate)

    ! pdaf
    use PDAFomi_obs_f, only:  PDAFomi_gather_obsstate, debug
    use PDAF, only: PDAF_get_obsmemberid

    ! intern
    use georeferenced_data_module, only: assign_to_observation_vector
    use grid_observation_module, only: field_bundle
    use quantity_computation_module, only: zg_mid_mean, zg_mid_mean_nm, zg_int_mean, zg_int_mean_nm
    use state_module, only:  state_vector

    IMPLICIT NONE

! *** Arguments ***
    class(obs_tme_grid) :: this
    INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
    INTEGER, INTENT(in) :: dim_obs               !< Dimension of full observed state (all observed fields)
    REAL, INTENT(in)    :: state_p(dim_p)        !< PE-local model state
    REAL, INTENT(inout) :: ostate(dim_obs)       !< Full observed state

! local
    REAL, ALLOCATABLE :: ostate_p(:)       ! local observed part of state vector
    type( field_bundle ) :: regridded
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
          ! TODO which ZG do we need to pass?
            call this%obs_op_computation(state_p, state_vector, regridded, zg_mid_mean, zg_mid_mean_nm, zg_int_mean, zg_int_mean_nm)
            call assign_to_observation_vector(this%tme_obs, regridded, ostate_p)
            call regridded%destroy

        else
          ! values have already bean interpolated in collect_state_pdaf.F90
          ostate_p = this%computed_observations(:,current_model_instance)
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

  END SUBROUTINE obs_tme_grid_obs_op

  !-------------------------------------------------------------------------------
!> Implementation of observation operator
!!
!! This routine applies the full observation operator
!! for the type of observations handled in this module.
!!
!! One can choose a proper observation operator from
!! PDAFOMI_OBS_OP or add one to that module or
!! implement another observation operator here.
!!
!! The routine is called by all filter processes.
!!
  subroutine obs_tme_grid_obs_op_computation(this, state_p, state_map, regridded, zg_mid, zg_mid_nm, zg_int, zg_int_nm)

    ! extern
    use ESMF, only: ESMF_Field

    ! tie-gcm
    use fields_module,only: levd0,levd1,lond0,lond1,latd0,latd1

    ! intern
    use configuration, only: cfg_filter
    use grid_observation_module, only: field_bundle, construct_field_bundle
    use quantity_info_module, only: quantity_info, get_info
    use state_module, only:  state_vector_mapping
    use tiegcm_optimized_interpolator, only: state_interpolator

    implicit none

    ! ATTENTION remember to deconstruct regridded after calling this function

    class(obs_tme_grid) :: this
    real, dimension(:), intent(in) :: state_p
    type(state_vector_mapping), intent(in) :: state_map
    type( field_bundle ), intent(inout) :: regridded
    real, dimension(levd0:levd1,lond0:lond1,latd0:latd1), intent(inout), optional :: zg_mid, zg_mid_nm, zg_int, zg_int_nm

    ! local
    type(state_interpolator) :: interp
    integer :: i
    type(quantity_info), pointer :: quantity
    type(ESMF_Field),pointer :: field

    regridded = construct_field_bundle(this%tme_obs%root%grid, this%field_names)

    ! 'If you pass an omitted dummy argument as the actual argument to a procedure, the corresponding dummy argument is considered to be omitted as well' (https://www.intel.com/content/www/us/en/develop/documentation/fortran-compiler-oneapi-dev-guide-and-reference/top/language-reference/program-units-and-procedures/argument-association-in-procedures/optional-arguments.html)
    call interp%init(this%tme_obs%bundle%bundle(1), &
                     zg_mid=zg_mid,&
                     zg_mid_nm=zg_mid_nm,&
                     zg_int=zg_int,&
                     zg_int_nm=zg_int_nm,&
                     degree=cfg_filter%spline_degree)

    ! loop over all fields
    do i = lbound(this%field_names,dim=1), ubound(this%field_names,dim=1)
      quantity=>get_info(this%field_names(i))
      call quantity%calc(state_map, state_p)
      field=>regridded%field(name=this%field_names(i))
      call interp%interpolate(quantity,quantity%data,field)
    end do

    call interp%deallocate()

  end subroutine obs_tme_grid_obs_op_computation

  !-------------------------------------------------------------------------------
!> Initialize local information on the module-type observation
!!
!! The routine is called during the loop over all local
!! analysis domains. It has to initialize the information
!! about local observations of the module type. It returns
!! number of local observations of the module type for the
!! current local analysis domain in DIM_OBS_L and the full
!! and local offsets of the observation in the overall
!! observation vector.
!!
!! This routine calls the routine PDAFomi_init_dim_obs_l
!! for each observation type. The call allows to specify a
!! different localization radius and localization functions
!! for each observation type and  local analysis domain.
!!
  SUBROUTINE obs_tme_grid_init_dim_obs_l(this, domain_p, step, dim_obs, dim_obs_l)

    ! PDAF
    use PDAF, only: PDAFomi_init_dim_obs_l
    use PDAFomi_obs_l, only: obs_l

    ! intern
    use mod_assimilation, only: coords_l, cutoff_radius, locweight, support_radius

    IMPLICIT NONE

! *** Arguments ***
    class(obs_tme_grid) :: this
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
    ! extent, since PDAF requires SIZE(cradius) == ncoord.

    CALL PDAFomi_init_dim_obs_l(local_obs, this%full_obs, &
         coords_l(1:this%full_obs%ncoord), locweight, &
         cutoff_radius(1:this%full_obs%ncoord), &
         support_radius(1:this%full_obs%ncoord), dim_obs_l)

  END SUBROUTINE obs_tme_grid_init_dim_obs_l



!-------------------------------------------------------------------------------
!> Perform covariance localization for local EnKF on the module-type observation
!!
!! The routine is called in the analysis step of the localized
!! EnKF. It has to apply localization to the two matrices
!! HP and HPH of the analysis step for the module-type
!! observation.
!!
!! This routine calls the routine PDAFomi_localize_covar
!! for each observation type. The call allows to specify a
!! different localization radius and localization functions
!! for each observation type.
!!
  SUBROUTINE obs_tme_grid_localize_covar(this, dim_p, dim_obs, HP_p, HPH, coords_p)

    ! PDAF
    use PDAF, only: PDAFomi_localize_covar

    ! intern
    use mod_assimilation, only: cutoff_radius, locweight, support_radius

    IMPLICIT NONE

! *** Arguments ***
    class(obs_tme_grid) :: this
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

  END SUBROUTINE obs_tme_grid_localize_covar

  !> Initializes the per-process observation coordinate array in spherical or TIE-GCM cell-index coordinates, as required by the configured localization filter.
  subroutine obs_tme_grid_init_coordinate_arrays(this, dim_obs_p, ocoord_p)

    ! intern
    use cell_id_coordinate_system, only: to_zonal_meridional_vertical
    use configuration, only: cfg_filter
    use mod_assimilation, only: filtertype
    use tgcm_pdaf_omi_obs_type_module, only: COORD_SPH_3D, COORD_CELL_IDX_3D
    use grid_observation_module, only: get_ocoord_lon_lat_lev, get_ocoord_tiegcm_grid_cell_idx_lin
    use quantity_info_module, only: LEVEL_MID, LEVEL_INT, LEVEL_NONE
    use quantity_computation_module, only: zg_mid_mean, zg_int_mean

    implicit none

    ! arguments
    class(obs_tme_grid) :: this

    INTEGER, INTENT(in) :: dim_obs_p    !< Dimension of full observation vector
    REAL, ALLOCATABLE, intent(out) :: ocoord_p(:,:)   ! PE-local observation coordinates

    ! local

    ! *** Initialize coordinate array of observations on the process sub-domain ***

    ALLOCATE(ocoord_p(this%full_obs%ncoord, dim_obs_p))

    if(filtertype==7) then

      select case(cfg_filter%localization_coord_sys)
        case(COORD_SPH_3D)
          ! quantities without vertical extent (LEVEL_NONE), i.e. VTEC, are
          ! localized horizontally only. ocoord_p then has two rows and
          ! get_ocoord_lon_lat_lev omits the altitude, see
          ! init_full_observation_dimension in obs_tgcm_pdaf_omi.F90.
          call get_ocoord_lon_lat_lev(ocoord_p,this%tme_obs)
        case(COORD_CELL_IDX_3D)
          ! ATTENTION only valid if all fields assimilated in this routine are
          ! located at mid or int points
          select case(this%tiegcm_native_grid)
            case(LEVEL_MID)
              call get_ocoord_tiegcm_grid_cell_idx_lin(zg_mid_mean,LEVEL_MID,ocoord_p,this%tme_obs)
              call to_zonal_meridional_vertical(ocoord_p)
            case(LEVEL_INT)
              call get_ocoord_tiegcm_grid_cell_idx_lin(zg_int_mean,LEVEL_INT,ocoord_p,this%tme_obs)
              call to_zonal_meridional_vertical(ocoord_p)
            case(LEVEL_NONE)
              ! TODO horizontal-only localization is implemented for spherical
              ! coordinates only. The cell index of a column-integrated
              ! quantity would need a lon/lat-only variant of
              ! get_ocoord_tiegcm_grid_cell_idx_lin, which derives the vertical
              ! cell index from a level.
              call shutdown('obs_tme_grid_init_coordinate_arrays: quantities without vertical'// &
                            ' extent can only be localized in spherical coordinates'// &
                            ' (localization_coord_sys=1). Observation: '//trim(this%name))
            case default
              ! without this the coordinates would silently stay uninitialized
              call shutdown('obs_tme_grid_init_coordinate_arrays: unhandled level of observation '//trim(this%name))
          end select
        case default
          call shutdown('invalid choice for localization coordinate system')
      end select

    else
      ocoord_p=0
      write(*,*) "ocoord_p are initialized with zeros", &
            " since they are not require for global filters"
    end if
  end subroutine

end module obs_tme_grid_pdafomi

MODULE obs_den_grid_pdafomi

  ! intern
  use configuration, only: cal_den_tme_file
  use grid_observation_module, only: reg_grid_dataset_group, reg_grid_dataset_root
  use obs_tme_grid_pdafomi, only: obs_tme_grid
  use pdaf_omi_obs_type_module, only: observation_interface

  IMPLICIT NONE

  type, extends(obs_tme_grid) :: obs_den_grid

    type(reg_grid_dataset_group) :: tme_distance

    type(cal_den_tme_file) :: config

    integer :: obs_out_ncid

    contains
      ! derived procedures
      procedure, pass(this) :: init_dim_obs => obs_den_grid_init_dim_obs
      procedure, pass(this) :: deallocate => obs_den_grid_deallocate
      procedure, pass(this) :: prodRinvA_l => obs_den_grid_prodRinvA_l
      ! Implementation specific procedures
      procedure, pass(this) :: init => obs_den_grid_init
      procedure, pass(this), private :: weights_from_distance_to_sat =>obs_den_grid_weights_from_distance_to_sat
  end type

!-------------------------------------------------------------------------------

CONTAINS

!> Opens the density observation, standard-deviation, and satellite-distance datasets for this density-grid observation type and links it to the result writer.
subroutine obs_den_grid_init(this, config)

  ! extern
  use netcdf

  ! tie-gcm
  use mpi_module, only: mytid

  ! intern
  use configuration, only: cfg_output
  use grid_observation_module, only: construct_reg_grid_dataset_group, construct_reg_grid_dataset_root
  use mod_parallel_pdaf, only: n_modeltasks, rank_world
  use netcdf_functionality, only: add_global_meta_data
  use quantity_info_module, only: LEVEL_MID

  implicit none

  ! arguments
  class(obs_den_grid) :: this
  type(cal_den_tme_file), intent(in) :: config

  ! local
  integer :: istat
  logical, dimension(10) :: mask

  this%tiegcm_native_grid = LEVEL_MID

  this%config=config
  this%name = config%name

  this%assimilate = config%apply

  mask = config%fields/=""
  allocate(this%field_names(count(mask)))
  this%field_names=pack(config%fields,mask)

  this%tme_data_file = construct_reg_grid_dataset_root(nc_file=config%tme_grid_file,&
                                                       altmin=config%lb_height,&
                                                       altmax=config%ub_height)

  this%tme_obs = construct_reg_grid_dataset_group(this%tme_data_file, &
                                        grp_name="nrlmsis2.0", &
                                        field_names=this%field_names, &
                                        grid_name=trim(config%name)//" observation")
  this%tme_std = construct_reg_grid_dataset_group(this%tme_data_file, &
                                        grp_name="nrlmsis2.0", &
                                        field_names=(/"DEN_STD"/), &
                                        grid_name=trim(config%name)//" standard deviation")
  this%tme_distance = construct_reg_grid_dataset_group(this%tme_data_file, &
                                      grp_name="Distances", &
                                      field_names=(/trim(this%config%satellite)//"_sph_dist ",&
                                                    trim(this%config%satellite)//"_vert_dist"/), &
                                      grid_name=trim(config%name)//" distance")


  call this%spatial_domain%init(domain_name=trim(config%name),&
                                save_members=cfg_output%save_members,&
                                kmax=cfg_output%max_moment,&
                                write_every_sec=config%write_every_sec,&
                                force_write_on_update=config%force_write_on_update,&
                                save_n_steps_after_update=cfg_output%save_n_steps_after_update)

  call this%spatial_domain%link_dataset(this%tme_data_file)

  call this%dst%init(this%tme_data_file,this%spatial_domain,this%name)

  call this%link_to_writer(this%spatial_domain, this%field_names)

  if ((this%assimilate) .and. (rank_world == 0)) then
    istat = nf90_create(path='results_'//trim(config%name)//'_observations.nc',&
                      cmode=NF90_NETCDF4,&
                      ncid=this%obs_out_ncid)
    call add_global_meta_data(this%obs_out_ncid)
  end if



  if ((this%assimilate) .or. (config%always_save)) then
    ! the interpolation is already done in collect state, so that all processes can do it instead of
    ! sequential execution here on filter prcocesses
    allocate( this%computed_observations(this%tme_obs%map%size_R(mytid), n_modeltasks ) )
  end if

end subroutine obs_den_grid_init

!> Destroys this density-grid observation's datasets and output file, and frees its arrays.
subroutine obs_den_grid_deallocate(this)

  ! extern
  use netcdf

  ! intern
  use mod_parallel_pdaf, only: rank_world

  implicit none

  class(obs_den_grid) :: this

  integer :: istat

  call this%tme_obs%destroy
  call this%tme_std%destroy
  call this%tme_distance%destroy
  call this%tme_data_file%destroy
  if (allocated (this%computed_observations))   deallocate (this%computed_observations)
  if (allocated (this%field_names))   deallocate (this%field_names)
  if ((this%assimilate) .and. (rank_world == 0)) then
    istat = nf90_close(this%obs_out_ncid)
  end if

  call this%dst%deallocate()

end subroutine obs_den_grid_deallocate

!> Initialize information on the module-type observation
!!
!! The routine is called by each filter process.
!! at the beginning of the analysis step before
!! the loop through all local analysis domains.
!!
!! It has to count the number of observations of the
!! observation type handled in this module according
!! to the current time step for all observations
!! required for the analyses in the loop over all local
!! analysis domains on the PE-local state domain.
!!
!! The following four variables have to be initialized in this routine
!! * full_obs\%doassim     - Whether to assimilate this type of observations
!! * full_obs\%disttype    - type of distance computation for localization with this observaton
!! * full_obs\%ncoord      - number of coordinates used for distance computation
!! * full_obs\%id_obs_p    - index of module-type observation in PE-local state vector
!!
!! Optional is the use of
!! * full_obs\%icoeff_p    - Interpolation coefficients for obs. operator (only if interpolation is used)
!! * full_obs\%domainsize  - Size of domain for periodicity for disttype=1 (<0 for no periodicity)
!! * full_obs\%obs_err_type - Type of observation errors for particle filter and NETF (default: 0=Gaussian)
!! * full_obs\%use_global obs - Whether to use global observations or restrict the observations to the relevant ones
!!                          (default: 1=use global full observations)
!!
!! Further variables are set when the routine PDAFomi_gather_obs is called.
!!
!! **Adapting the template**
!! In this routine the variables listed above have to be initialized. One
!! can include modules from the model with 'use', e.g. for mesh information.
!! Alternatively one could include these as subroutine arguments
!!
  SUBROUTINE obs_den_grid_init_dim_obs(this, step, dim_obs)

    ! extern
    use esmf

    ! tie-gcm
    use mpi_module, only: mytid
    use params_module, only : nlon

    ! pdaf
    use PDAF, only: PDAFomi_gather_obs

    ! intern
    use array_print_module, only: printMat
    use array_mapping_module, only: flatten, map_to_3d
    use mod_assimilation, only: cutoff_radius, analysis_step_count
    use result_file_writer_module, only: instant_reg_grid_write

    IMPLICIT NONE

! *** Arguments ***
    class(obs_den_grid) :: this
    INTEGER, INTENT(in)    :: step       !< Current time step
    INTEGER, INTENT(inout) :: dim_obs    !< Dimension of full observation vector

! *** Local variables ***
    INTEGER :: i
    INTEGER :: dim_obs_p                 ! Number of process-local observations
    REAL, ALLOCATABLE :: obs_p(:)        ! PE-local observation vector
    REAL, ALLOCATABLE :: ivar_obs_p(:)   ! PE-local inverse observation error variance
    REAL, ALLOCATABLE :: ocoord_p(:,:)   ! PE-local observation coordinates

    real(ESMF_KIND_R8), pointer :: field_data(:,:,:)
    integer :: rc

    real, allocatable, dimension(:,:,:) :: dist_weights, std_den
    real, dimension(:,:,:), contiguous, pointer :: mat3d

!     real, allocatable ::debug(:,:)


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************
  call this%init_full_observation_dimension

! **********************************
! *** Read PE-local observations ***
! **********************************

  ! read observation values and their coordinates
  ! also read observation error information if available
  call this%tme_obs%read_current_epoch
  call this%tme_std%read_current_epoch
  call this%tme_distance%read_current_epoch
! ***********************************************************
! *** Count available observations for the process domain ***
! *** and initialize index and coordinate arrays.         ***
! ***********************************************************

    ! *** Count valid observations that lie within the process sub-domain ***

    dim_obs_p = this%tme_obs%map%size_R(mytid)

    ! *** Initialize vector of observations on the process sub-domain ***

    ALLOCATE(obs_p(dim_obs_p))

    do i = 1, this%tme_obs%bundle%size()
      call ESMF_FieldGet(this%tme_obs%bundle%field(idx=i), localDe=0, farrayPtr=field_data, &
            rc=rc)
      obs_p( this%tme_obs%map%idx_R(i,mytid)%begin_p : &
             this%tme_obs%map%idx_R(i,mytid)%back_p ) &
           = RESHAPE( field_data, (/this%tme_obs%map%F_R_size(i,mytid)/) )

    end do

    ! *** Initialize coordinate array of observations on the process sub-domain ***
    call this%init_coordinate_arrays(dim_obs_p, ocoord_p)

    ! *** Initialize process local index array                         ***
    ! *** This array holds the information which elements of the state ***
    ! *** vector are used in the observation operator.                 ***
    ! *** It has a many rows as required for the observation operator, ***
    ! *** i.e. 1 if observations are at grid points; >1 if             ***
    ! *** interpolation is required                                    ***

  ! The initialization is done locally for each process sub-domain and later
  ! used in the observation operator.
  ! Examples:
  ! 1. If the observations are model fields located at grid points, one should
  !   initialize the index array full_obs%id_obs_p with one row so that it contains
  !   the indices of the observed field values in the process-local state vector
  !   (state_p). Then one can use the observation operator OBS_OP_GRIDPOINT
  !   provided by the module PDAF.
  ! 2. If the observations are the average of model fields located at grid points,
  !   one should initialize the index array full_obs%id_obs_p with as many rows as
  !   values to be averaged. Each column of the arrays then contains the indices of
  !   the elements of the process-local state vector that have to be averaged. With
  !   this index array one can then use the observation operator OBS_OP_GRIDAVG
  !   provided by the module PDAF.
  ! 3. If model values need to be interpolated to the observation location
  !   one should initialize the index array full_obs%id_obs_p with as many rows as
  !   values are required in the interpolationto be averaged. Each column of the
  !   array then contains the indices of elements of the process-local state vector
  !   that are used in the interpolation.
  ! Below, you need to replace NROWS by the number of required rows

!    ALLOCATE(full_obs%id_obs_p( NROWS , dim_obs_p))
    ALLOCATE(this%full_obs%id_obs_p(0,0))
!    full_obs%id_obs_p = ...


! **********************************************************************
! *** Initialize interpolation coefficients for observation operator ***
! **********************************************************************

  ! This initialization is only required if an observation operator
  ! with interpolation is used. The coefficients should be determined
  ! here instead of the observation operator, because the operator is
  ! called for each ensemble member while init_dim_obs is only called
  ! once.

  ! Allocate array of interpolation coefficients. As full_obs%id_obs_p, the number
  ! of rows corresponds to the number of grid points using the the interpolation

!    ALLOCATE(full_obs%icoeff_p( NROWS , dim_obs_p))

  ! Ensure that the order of the coefficients is consistent with the
  ! indexing in full_obs%id_obs_p. Further ensure that the order is consistent
  ! with the assumptions used in the observation operator.

!    full_obs%icoeff_p = ...

! ****************************************************************
! *** Define observation errors for process-local observations ***
! ****************************************************************

    ALLOCATE(ivar_obs_p(dim_obs_p))

    do i = 1, this%tme_obs%bundle%size()
      if((this%tme_obs%bundle%names(i) == "DEN") .or. (this%tme_obs%bundle%names(i) == "DEN_NM")  ) then

        allocate( std_den(this%tme_obs%root%nlon,this%tme_obs%root%nlat,this%tme_obs%root%nalt))
        call this%weights_from_distance_to_sat(dist_weights)

        call ESMF_FieldGet(this%tme_std%bundle%field(name="DEN_STD"), localDe=0, farrayPtr=field_data, &
              rc=rc)

        std_den(:,:,:) = field_data/(dist_weights*this%config%weight)

        ivar_obs_p( this%tme_obs%map%idx_R(i,mytid)%begin_p : &
                    this%tme_obs%map%idx_R(i,mytid)%back_p ) &
                  = RESHAPE( std_den, (/this%tme_obs%map%F_R_size(i,mytid)/) )

        call instant_reg_grid_write(this%obs_out_ncid,"dist_weights",analysis_step_count,dist_weights,this%tme_data_file)

        deallocate(dist_weights,std_den)
      else
        ivar_obs_p( this%tme_obs%map%idx_R(i,mytid)%begin_p : &
                    this%tme_obs%map%idx_R(i,mytid)%back_p ) &
           = obs_p( this%tme_obs%map%idx_R(i,mytid)%begin_p : &
                    this%tme_obs%map%idx_R(i,mytid)%back_p )/10
      end if

      ! write std
      call map_to_3d(flat=ivar_obs_p,&
                      mat3d=mat3d,&
                      mat_ub=(/this%tme_obs%root%nlon,this%tme_obs%root%nlat,this%tme_obs%root%nalt/),&
                      first=this%tme_obs%map%idx_R(i,mytid)%begin_p,&
                      last=this%tme_obs%map%idx_R(i,mytid)%back_p)

      call instant_reg_grid_write(this%obs_out_ncid,&
                                  'std_'//this%tme_obs%bundle%names(i),&
                                  analysis_step_count,&
                                  mat3d,&
                                  this%tme_data_file)
      ! write obs
      call map_to_3d(flat=obs_p,&
                      mat3d=mat3d,&
                      mat_ub=(/this%tme_obs%root%nlon,this%tme_obs%root%nlat,this%tme_obs%root%nalt/),&
                      first=this%tme_obs%map%idx_R(i,mytid)%begin_p,&
                      last=this%tme_obs%map%idx_R(i,mytid)%back_p)

      call instant_reg_grid_write(this%obs_out_ncid,&
                                  this%tme_obs%bundle%names(i),&
                                  analysis_step_count,&
                                  mat3d,&
                                  this%tme_data_file)

    end do

    ivar_obs_p = 1.0/ivar_obs_p**2

! ****************************************
! *** Gather global observation arrays ***
! ****************************************

    ! NOTE FOR DIM_OBS_P=0
    ! For the call to PDAFomi_gather_obs_f, obs_p, ivar_obs_p, ocoord_p,
    ! and full_obs%id_obs_p need to be allocated. Thus, if dim_obs_p=0 can
    ! happen in your application you should explicitly handle this case.
    ! You can introduce an IF block in the initializations above:
    !  IF dim_obs_p>0 THEN
    !     regular allocation and initialization of obs_p, ivar_obs_p, ocoord_p
    !  ELSE
    !     allocate obs_p, ivar_obs_p, ocoord_p, full_obs%id_obs_p with size=1
    !  ENDIF


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

! ***************
! TEST OUTPUT
! ***************
!     allocate(debug(dim_obs_p,5))
!     debug(:,1:3) = transpose(ocoord_p)
!     debug(:,4) = obs_p
!     debug(:,5) = sqrt(1/ivar_obs_p)
!
!     call printMat(debug, "observation table")
!     deallocate(debug)


! ********************
! *** Finishing up ***
! ********************

    ! Deallocate all local arrays
    DEALLOCATE(obs_p, ocoord_p, ivar_obs_p)

    ! Arrays in THISOBS have to be deallocated after the analysis step
    ! by a call to deallocate_obs() in prepoststep_pdaf.

  END SUBROUTINE obs_den_grid_init_dim_obs

  !> Computes horizontal/vertical exponential-decay weights from each grid point's distance to the satellite track.
  SUBROUTINE obs_den_grid_weights_from_distance_to_sat(this, P)

    ! extern
    use esmf, only: ESMF_FieldGet, ESMF_KIND_R8

    ! intern
    use array_print_module, only: printMat

    IMPLICIT NONE

      ! arguments
      class(obs_den_grid) :: this
      real, allocatable, dimension(:,:,:), intent(inout) :: P

      ! local
      real(ESMF_KIND_R8), dimension(:,:,:), pointer :: sph_dist, vert_dist

      !
      ! maximal spherical distance is 180 degrees:
      !   half life  | weight(45) | weight(90) | weight(180)
      !  180   deg   |     0.84   |  1/sqrt(2) |    1/2
      !   90   deg   |  1/sqrt(2) |     1/2    |    1/4
      !   45   deg   |     1/2    |     1/4    |    1/16
      !   22.5 deg   |     1/4    |     1/16   |    1/256

      real, parameter :: ln2 = log(2.)

      call ESMF_FieldGet(this%tme_distance%bundle%field(name=trim(this%config%satellite)//"_sph_dist"),&
                         localDe=0,&
                         farrayPtr=sph_dist)

      call ESMF_FieldGet(this%tme_distance%bundle%field(name=trim(this%config%satellite)//"_vert_dist"),&
                         localDe=0,&
                         farrayPtr=vert_dist)

      allocate(P, mold=sph_dist)

      ! exponential decay: f(x) = A0 exp(-lambda * x)
      !
      ! half life:  hw = ln(2)/lambda
      !
      ! P =  exp(-ln(2)/hw * x)
      !


      P = 1
      if(this%config%horz_weight_half_life .gt. 0) then
        P(:,:,:) = exp( -ln2/this%config%horz_weight_half_life * sph_dist)
      end if
      if(this%config%vert_weight_half_life .gt. 0) then
        P(:,:,:) = P*exp( -ln2/this%config%vert_weight_half_life * abs(vert_dist))
      end if

  END SUBROUTINE obs_den_grid_weights_from_distance_to_sat

  !> Computes C = R^-1 A for the local density observations, using either a diagonal or a full correlated observation-error covariance matrix.
  SUBROUTINE obs_den_grid_prodRinvA_l(this,domain_p, step, dim_obs_l, rank, obs_l, A_l, C_l)

    ! PDAF
    use PDAF, only: PDAFomi_prodRinvA_l, PDAFomi_observation_localization_weights
    ! the dummy argument obs_l below shadows the type name, hence the rename
    use PDAFomi_obs_l, only: obs_l_type => obs_l
    use default_lapack_interface, only: default_dposv

    ! intern
    use array_print_module,&
      only: printMat

    IMPLICIT NONE

    ! *** Arguments ***
    class(obs_den_grid) :: this
    INTEGER, INTENT(in) :: domain_p             !< Index of current local analysis domain
    INTEGER, INTENT(in) :: step                 !< Current time step
    INTEGER, INTENT(in) :: dim_obs_l            !< Dimension of local observation vector
    INTEGER, INTENT(in) :: rank                 !< Rank of initial covariance matrix
    REAL, INTENT(in)    :: obs_l(dim_obs_l)     !< Local vector of observations
    REAL, INTENT(inout) :: A_l(dim_obs_l, rank) !< Input matrix
    REAL, INTENT(out)   :: C_l(dim_obs_l, rank) !< Output matrix

    ! local
    type(obs_l_type), pointer :: local_obs ! this thread's local observation
    REAL, ALLOCATABLE :: weight(:)     ! Localization weights
    REAL, ALLOCATABLE :: R(:,:)
    integer :: i, j
    integer :: off

    logical, save :: first_call = .true.

    if(dim_obs_l<1) return

    call this%local_obs(local_obs)

    if(trim(this%config%correlations)=="none") then
      if(first_call) write(*,*) "use diagonal vcm"
      CALL PDAFomi_prodRinvA_l(local_obs, this%full_obs, dim_obs_l, rank, A_l, C_l, 0)
    else

      if(first_call) write(*,*) "use full vcm"

      off = local_obs%off_obs_l

      allocate(R(dim_obs_l,dim_obs_l))

      ! TODO specific for each omi type
      call calc_obs_l_vcm(local_obs, this%full_obs, R)

      if(local_obs%locweight/=0) then ! observation localization

        allocate(weight(dim_obs_l))

        call PDAFomi_observation_localization_weights(local_obs, this%full_obs, rank, A_l, &
                                                      weight, 0)

        ! (1) Weights must be non negative since we compute the sqrt.
        !     Exp and Gaspari and Cohn (GC) might give small negative results
        !     close to the end of the finite support due to numeric issues
        ! (2) We divide the covariance through the weight, thus the weights
        !     are not allowed to be zero
        where(weight<1E-15)
          weight=1E-15
        elsewhere
          weight = sqrt(weight)
        end where

        ! (W-^1 R W^-1)^-1 = W R^-1 W
        ! W is a diagonal matrix containing the square roots of the weights
        ! apply weights
        do j = 1, dim_obs_l
          do i = j, dim_obs_l
            R(i,j) = R(i,j)/(weight(i)*weight(j))
          end do
        end do

      end if

!       call printMat(A_l,'A_l')

      !!!!!!!!!!!!!!!!!!
      ! solve R C = A  !
      ! C = R^(-1) A   !
      !!!!!!!!!!!!!!!!!!
      C_l(off+1:off+dim_obs_l,:) = A_l(off+1:off+dim_obs_l,:)
      call default_dposv("L",R,C_l(off+1:off+dim_obs_l,:))

!       call printMat(C_l,'C_l')

      if(allocated(weight))deallocate(weight)
      deallocate(R)
    end if

    if(first_call) first_call = .false.

  END SUBROUTINE obs_den_grid_prodRinvA_l

  !> Builds the local observation error variance-covariance matrix (lower triangle) from observation standard deviations and a distance-based correlation function.
  subroutine calc_obs_l_vcm(thisobs_l, thisobs, R)

    ! intern
    use array_print_module,&
      only: printMat

    ! pdaf
    use PDAFomi_obs_l, only: obs_l
    use PDAFomi_obs_f, only: obs_f

    implicit none

    ! arguments
    TYPE(obs_l), INTENT(inout) :: thisobs_l  !< Data type with local observation
    TYPE(obs_f), INTENT(inout) :: thisobs    !< Data type with full observation
    real, dimension(thisobs_l%dim_obs_l,thisobs_l%dim_obs_l), intent(out) :: R

    ! local
    integer :: i
    integer :: off

    real, dimension(thisobs_l%dim_obs_l) :: sigma
    real, dimension(thisobs_l%dim_obs_l,thisobs%ncoord) :: coordinates

    off = thisobs_l%off_obs_l

    ! calcuate standard deviation
    do i=1,size(R,dim=1)
      sigma(i) = 1./sqrt(thisobs_l%ivar_obs_l(off+i))
    end do

    do i=1, size(R,dim=1)
        coordinates(i,:) =  thisobs%ocoord_f(:, thisobs_l%id_obs_l(off+i))
    end do

    call calc_correlation_matrix(thisobs_l%dim_obs_l,coordinates,R)

    call vcm_from_std_and_corr(thisobs_l%dim_obs_l,sigma,R)

  end subroutine calc_obs_l_vcm

  !> Fills the lower triangle of a correlation matrix from an exponential decay in great-circle and vertical distance between observation points.
  subroutine calc_correlation_matrix(n,coordinates,C)

    implicit none

    ! arguments
    integer, intent(in) :: n
    real, dimension(n,3), intent(in) :: coordinates
    real, dimension(n,n), intent(out) :: C

    ! local
    integer :: i, j

    real :: great_circle_dist
    real :: vertical_dist
    real :: slon, slat

    real :: dist

    C = 0
    ! fill main diagonal with ones
    do i=1,size(C,dim=1)
      C(i,i) = 1
    end do

    ! fill lower triangle with correlations
    do j=1, size(C,dim=2)
      do i=j+1, size(C,dim=1)


        slon = SIN((coordinates(i,1) - coordinates(j,1))/2)
        slat = SIN((coordinates(i,2) - coordinates(j,2))/2)

        great_circle_dist = SQRT(slat*slat + COS(coordinates(i,2))*COS(coordinates(j,2))*slon*slon)
        great_circle_dist = 2.0 * 6371 * ASIN(great_circle_dist) ! [km]
        vertical_dist = ABS(coordinates(i,3) - coordinates(j,3))/1000 ! [km]

        dist = sqrt((great_circle_dist/2670)**2+vertical_dist**2)
        C(i,j) = exp(-dist/20.5)

!       write(*,*) vertical_dist, great_circle_dist, C(i,j)

      end do
    end do

  !    call printMat(R,'local obs correlation matrix')

  end subroutine

 !> Converts a correlation matrix in-place to a variance-covariance matrix (VCM = diag(S) R diag(S), lower triangle only).
 !! TODO move this subroutine to a more general location
 subroutine vcm_from_std_and_corr(n,sigma,R)

    ! intern
    use array_print_module,&
      only: printMat

    implicit none

    ! arguments
    integer, intent(in) :: n
    real, dimension(n), intent(in) :: sigma
    real, dimension(n,n), intent(inout) :: R ! on entry correlation in lower trinagle, on exit vcm in lower triangle

    ! local
    integer :: i, j

    if(n<1)then ! no observation
      return
    elseif(n==1)then ! one observation
      ! R is 1x1 and contains only the variance of the observation
      R(1,1) = sigma(1)*sigma(1)
    else
      do j=1, size(R,dim=2)
        do i=j, size(R,dim=1)
          R(i,j) = R(i,j) * sigma(i) * sigma(j)
        end do
      end do
    end if

!     call printMat(R,'local obs vcm')

  end subroutine vcm_from_std_and_corr

END MODULE obs_den_grid_pdafomi



MODULE obs_tum_ne_pdafomi

  ! intern
  use configuration, only: tum_ne_tme_file
  use grid_observation_module, only: reg_grid_dataset_group, reg_grid_dataset_root
  use obs_tme_grid_pdafomi, only: obs_tme_grid
  use pdaf_omi_obs_type_module, only: observation_interface
  use result_file_writer_module, only: nc_reg_grid

  IMPLICIT NONE

  type, extends(obs_tme_grid) :: obs_tum_ne
    type(tum_ne_tme_file) :: config
    integer :: obs_out_ncid
    contains
      ! derived procedures
      procedure, pass(this) :: init_dim_obs => obs_tum_ne_init_dim_obs
      procedure, pass(this) :: deallocate => obs_tum_ne_deallocate
      ! Implementation specific procedures
      procedure, pass(this) :: init => obs_tum_ne_init
  end type

!-------------------------------------------------------------------------------

CONTAINS

!> Opens the TUM electron-density observation and standard-deviation datasets and links this observation type to the result writer.
subroutine obs_tum_ne_init(this, config)

  ! extern
  use netcdf

  ! tie-gcm
  use mpi_module, only: mytid

  ! intern
  use configuration, only: cfg_output
  use grid_observation_module, only: construct_reg_grid_dataset_group, construct_reg_grid_dataset_root
  use mod_parallel_pdaf, only: n_modeltasks, rank_world
  use netcdf_functionality, only: add_global_meta_data
  use quantity_info_module, only: LEVEL_INT

  implicit none

  class(obs_tum_ne) :: this
  type(tum_ne_tme_file), intent(in) :: config

  integer :: istat

  this%tiegcm_native_grid = LEVEL_INT

  this%config=config
  this%name = config%name

  this%assimilate = config%apply

  allocate(this%field_names(1))
  this%field_names = (/"NE"/)

  this%tme_data_file = construct_reg_grid_dataset_root(nc_file=config%tme_grid_file,altmin=config%lb_height,altmax=config%ub_height)

  this%tme_obs = construct_reg_grid_dataset_group(this%tme_data_file, &
                                        field_names=this%field_names, &
                                        grid_name=trim(config%name)//" observation")
  this%tme_std = construct_reg_grid_dataset_group(this%tme_data_file, &
                                        field_names=(/"NE_STD"/), &
                                        grid_name=trim(config%name)//" standard deviation")

  call this%spatial_domain%init(domain_name=trim(config%name),&
                                save_members=cfg_output%save_members,&
                                kmax=cfg_output%max_moment,&
                                write_every_sec=config%write_every_sec,&
                                force_write_on_update=config%force_write_on_update,&
                                save_n_steps_after_update=cfg_output%save_n_steps_after_update)

  call this%spatial_domain%link_dataset(this%tme_data_file)
  call this%dst%init(this%tme_data_file,this%spatial_domain,this%name)
  call this%link_to_writer(this%spatial_domain, this%field_names)

  if (rank_world == 0) then
    istat = nf90_create(path='results_'//trim(config%name)//'_observations.nc',&
                      cmode=NF90_NETCDF4,&
                      ncid=this%obs_out_ncid)
    call add_global_meta_data(this%obs_out_ncid)
  end if

  if ((this%assimilate) .or. (config%always_save)) then
    ! the interpolation is already done in collect state, so that all processes can do it instead of
    ! sequential execution here on filter prcocesses
    allocate( this%computed_observations(this%tme_obs%map%size_R(mytid), n_modeltasks ) )
  end if

end subroutine obs_tum_ne_init

!> Destroys this TUM electron-density observation's datasets and output file, and frees its arrays.
subroutine obs_tum_ne_deallocate(this)

  ! extern
  use netcdf

  ! intern
  use mod_parallel_pdaf, only: rank_world

  implicit none

  class(obs_tum_ne) :: this

  integer :: istat

  call this%tme_obs%destroy
  call this%tme_std%destroy
  call this%tme_data_file%destroy
  if (allocated (this%computed_observations))   deallocate (this%computed_observations)
  if (allocated (this%field_names))   deallocate (this%field_names)
  if (rank_world == 0) then
    istat = nf90_close(this%obs_out_ncid)
  end if

  call this%dst%deallocate()

end subroutine obs_tum_ne_deallocate

!> Initialize information on the module-type observation
!!
!! The routine is called by each filter process.
!! at the beginning of the analysis step before
!! the loop through all local analysis domains.
!!
!! It has to count the number of observations of the
!! observation type handled in this module according
!! to the current time step for all observations
!! required for the analyses in the loop over all local
!! analysis domains on the PE-local state domain.
!!
!! The following four variables have to be initialized in this routine
!! * full_obs\%doassim     - Whether to assimilate this type of observations
!! * full_obs\%disttype    - type of distance computation for localization with this observaton
!! * full_obs\%ncoord      - number of coordinates used for distance computation
!! * full_obs\%id_obs_p    - index of module-type observation in PE-local state vector
!!
!! Optional is the use of
!! * full_obs\%icoeff_p    - Interpolation coefficients for obs. operator (only if interpolation is used)
!! * full_obs\%domainsize  - Size of domain for periodicity for disttype=1 (<0 for no periodicity)
!! * full_obs\%obs_err_type - Type of observation errors for particle filter and NETF (default: 0=Gaussian)
!! * full_obs\%use_global obs - Whether to use global observations or restrict the observations to the relevant ones
!!                          (default: 1=use global full observations)
!!
!! Further variables are set when the routine PDAFomi_gather_obs is called.
!!
!! **Adapting the template**
!! In this routine the variables listed above have to be initialized. One
!! can include modules from the model with 'use', e.g. for mesh information.
!! Alternatively one could include these as subroutine arguments
!!
  SUBROUTINE obs_tum_ne_init_dim_obs(this, step, dim_obs)

    ! extern
    use esmf

    ! tie-gcm
    use mpi_module, only: mytid
    use params_module, only : nlon

    ! pdaf
    use PDAF, only: PDAFomi_gather_obs

    ! intern
    use array_print_module, only: printMat
    use array_mapping_module, only: flatten, map_to_3d
!     use mod_assimilation, only: filtertype
    use mod_assimilation, only: cutoff_radius, analysis_step_count
    use result_file_writer_module, only: instant_reg_grid_write

    IMPLICIT NONE

! *** Arguments ***
    class(obs_tum_ne) :: this
    INTEGER, INTENT(in)    :: step       !< Current time step
    INTEGER, INTENT(inout) :: dim_obs    !< Dimension of full observation vector

! *** Local variables ***
    INTEGER :: i
    INTEGER :: dim_obs_p                 ! Number of process-local observations
    REAL, ALLOCATABLE :: obs_p(:)        ! PE-local observation vector
    REAL, ALLOCATABLE :: ivar_obs_p(:)   ! PE-local inverse observation error variance
    REAL, ALLOCATABLE :: ocoord_p(:,:)   ! PE-local observation coordinates

    real(ESMF_KIND_R8), pointer :: field_data(:,:,:)
    integer :: rc

    real, dimension(:,:,:), contiguous, pointer :: mat3d

!     real, allocatable ::debug(:,:)


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

  call this%init_full_observation_dimension


! **********************************
! *** Read PE-local observations ***
! **********************************

  ! read observation values and their coordinates
  ! also read observation error information if available
  call this%tme_obs%read_current_epoch
  call this%tme_std%read_current_epoch

! ***********************************************************
! *** Count available observations for the process domain ***
! *** and initialize index and coordinate arrays.         ***
! ***********************************************************

  ! *** Count valid observations that lie within the process sub-domain ***

  dim_obs_p = this%tme_obs%map%size_R(mytid)

  ! *** Initialize vector of observations on the process sub-domain ***

  ALLOCATE(obs_p(dim_obs_p))

  call ESMF_FieldGet(this%tme_obs%bundle%field(name='NE'), localDe=0, farrayPtr=field_data, &
          rc=rc)
  i=1
  obs_p( this%tme_obs%map%idx_R(i,mytid)%begin_p : &
         this%tme_obs%map%idx_R(i,mytid)%back_p ) &
        = RESHAPE( field_data, (/this%tme_obs%map%F_R_size(i,mytid)/) )

!   call printMat(field_data,name='obs',order=(/3,1,2/))

  ! convert 1/m3 to 1/cm3
  obs_p = obs_p/1E+6

  ! *** Initialize coordinate array of observations on the process sub-domain ***
  call this%init_coordinate_arrays(dim_obs_p, ocoord_p)

    ! *** Initialize process local index array                         ***
    ! *** This array holds the information which elements of the state ***
    ! *** vector are used in the observation operator.                 ***
    ! *** It has a many rows as required for the observation operator, ***
    ! *** i.e. 1 if observations are at grid points; >1 if             ***
    ! *** interpolation is required                                    ***

  ! The initialization is done locally for each process sub-domain and later
  ! used in the observation operator.
  ! Examples:
  ! 1. If the observations are model fields located at grid points, one should
  !   initialize the index array full_obs%id_obs_p with one row so that it contains
  !   the indices of the observed field values in the process-local state vector
  !   (state_p). Then one can use the observation operator OBS_OP_GRIDPOINT
  !   provided by the module PDAF.
  ! 2. If the observations are the average of model fields located at grid points,
  !   one should initialize the index array full_obs%id_obs_p with as many rows as
  !   values to be averaged. Each column of the arrays then contains the indices of
  !   the elements of the process-local state vector that have to be averaged. With
  !   this index array one can then use the observation operator OBS_OP_GRIDAVG
  !   provided by the module PDAF.
  ! 3. If model values need to be interpolated to the observation location
  !   one should initialize the index array full_obs%id_obs_p with as many rows as
  !   values are required in the interpolationto be averaged. Each column of the
  !   array then contains the indices of elements of the process-local state vector
  !   that are used in the interpolation.
  ! Below, you need to replace NROWS by the number of required rows

!    ALLOCATE(full_obs%id_obs_p( NROWS , dim_obs_p))
    ALLOCATE(this%full_obs%id_obs_p(0,0))

!    full_obs%id_obs_p = ...


! **********************************************************************
! *** Initialize interpolation coefficients for observation operator ***
! **********************************************************************

  ! This initialization is only required if an observation operator
  ! with interpolation is used. The coefficients should be determined
  ! here instead of the observation operator, because the operator is
  ! called for each ensemble member while init_dim_obs is only called
  ! once.

  ! Allocate array of interpolation coefficients. As full_obs%id_obs_p, the number
  ! of rows corresponds to the number of grid points using the the interpolation

!    ALLOCATE(full_obs%icoeff_p( NROWS , dim_obs_p))

  ! Ensure that the order of the coefficients is consistent with the
  ! indexing in full_obs%id_obs_p. Further ensure that the order is consistent
  ! with the assumptions used in the observation operator.

!    full_obs%icoeff_p = ...

! ****************************************************************
! *** Define observation errors for process-local observations ***
! ****************************************************************

      ALLOCATE(ivar_obs_p(dim_obs_p))

      call ESMF_FieldGet(this%tme_std%bundle%field(name="NE_STD"), localDe=0, farrayPtr=field_data, &
              rc=rc)

      ivar_obs_p(:) = RESHAPE( field_data, (/size(field_data)/) );

      ! convert 1/m3 to 1/cm3
      ivar_obs_p = ivar_obs_p/1E+6

      ! apply weight
      ivar_obs_p = ivar_obs_p/this%config%weight

      ! limit standard deviation for numeric stability
      ! tum ne (especially on the on night side) is very small
      ! observations is given to much trust
      where(ivar_obs_p<this%config%min_std) ivar_obs_p=this%config%min_std

      where (isnan(obs_p))
        ! nans occur in vincity of low electron density. set them to minimal density and downweight
        obs_p = 3100
        ivar_obs_p = 1E+9
      elsewhere (obs_p<3100 )
        ! limit smallest value to smallest value of TIE-GCM
        ! also down weight them,
        obs_p = 3100
        ivar_obs_p = 1E+7
      end where



      ! write std
      call map_to_3d(flat=ivar_obs_p,&
                      mat3d=mat3d,&
                      mat_ub=(/this%tme_obs%root%nlon,this%tme_obs%root%nlat,this%tme_obs%root%nalt/),&
                      first=this%tme_obs%map%idx_R(i,mytid)%begin_p,&
                      last=this%tme_obs%map%idx_R(i,mytid)%back_p)

     call instant_reg_grid_write(this%obs_out_ncid,&
                                  'std_'//this%tme_obs%bundle%names(i),&
                                  analysis_step_count,&
                                  mat3d,&
                                  this%tme_data_file)
      ! write obs
      call map_to_3d(flat=obs_p,&
                      mat3d=mat3d,&
                      mat_ub=(/this%tme_obs%root%nlon,this%tme_obs%root%nlat,this%tme_obs%root%nalt/),&
                      first=this%tme_obs%map%idx_R(i,mytid)%begin_p,&
                      last=this%tme_obs%map%idx_R(i,mytid)%back_p)

      call instant_reg_grid_write(this%obs_out_ncid,&
                                  this%tme_obs%bundle%names(i),&
                                  analysis_step_count,&
                                  mat3d,&
                                  this%tme_data_file)

     ! std to inverse variance
     ivar_obs_p = 1.0/(ivar_obs_p**2);

! ****************************************
! *** Gather global observation arrays ***
! ****************************************

    ! NOTE FOR DIM_OBS_P=0
    ! For the call to PDAFomi_gather_obs_f, obs_p, ivar_obs_p, ocoord_p,
    ! and full_obs%id_obs_p need to be allocated. Thus, if dim_obs_p=0 can
    ! happen in your application you should explicitly handle this case.
    ! You can introduce an IF block in the initializations above:
    !  IF dim_obs_p>0 THEN
    !     regular allocation and initialization of obs_p, ivar_obs_p, ocoord_p
    !  ELSE
    !     allocate obs_p, ivar_obs_p, ocoord_p, full_obs%id_obs_p with size=1
    !  ENDIF


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

! ***************
! TEST OUTPUT
! ***************
!     allocate(debug(dim_obs_p,5))
!     debug(:,1:3) = transpose(ocoord_p)
!     debug(:,4) = obs_p
!     debug(:,5) = sqrt(1/ivar_obs_p)
!
!     call printMat(debug, "observation table")
!     deallocate(debug)


! ********************
! *** Finishing up ***
! ********************

    ! Deallocate all local arrays
    DEALLOCATE(obs_p, ocoord_p, ivar_obs_p)

    ! Arrays in THISOBS have to be deallocated after the analysis step
    ! by a call to deallocate_obs() in prepoststep_pdaf.

  END SUBROUTINE obs_tum_ne_init_dim_obs

END MODULE obs_tum_ne_pdafomi


MODULE obs_tum_vtec_pdafomi

  ! intern
  use configuration, only: tum_vtec_file
  use grid_observation_module, only: reg_grid_dataset_group, reg_grid_dataset_root
  use obs_tme_grid_pdafomi, only: obs_tme_grid
  use pdaf_omi_obs_type_module, only: observation_interface
  use result_file_writer_module, only: nc_reg_grid

  IMPLICIT NONE

  type, extends(obs_tme_grid) :: obs_tum_vtec
    type(tum_vtec_file) :: config
    integer :: obs_out_ncid
    contains
      ! derived procedures
      procedure, pass(this) :: init_dim_obs => obs_tum_vtec_init_dim_obs
      procedure, pass(this) :: deallocate => obs_tum_vtec_deallocate
      ! VTEC has no vertical extent. Thus the observation operator of
      ! obs_tme_grid, which interpolates vertically against the geometric
      ! height, cannot be used and is overridden here.
      procedure, pass(this) :: obs_op_computation => obs_tum_vtec_obs_op_computation
      ! Implementation specific procedures
      procedure, pass(this) :: init => obs_tum_vtec_init
  end type

!-------------------------------------------------------------------------------

CONTAINS

!> Opens the TUM electron-density observation and standard-deviation datasets and links this observation type to the result writer.
subroutine obs_tum_vtec_init(this, config)

  ! extern
  use netcdf

  ! tie-gcm
  use mpi_module, only: mytid

  ! intern
  use configuration, only: cfg_output
  use grid_observation_module, only: construct_reg_grid_dataset_group, construct_reg_grid_dataset_root
  use mod_parallel_pdaf, only: n_modeltasks, rank_world
  use netcdf_functionality, only: add_global_meta_data
  use quantity_info_module, only: LEVEL_NONE

  implicit none

  class(obs_tum_vtec) :: this
  type(tum_vtec_file), intent(in) :: config

  integer :: istat

  this%tiegcm_native_grid = LEVEL_NONE
  this%config=config
  this%name = config%name

  this%assimilate = config%apply

  allocate(this%field_names(1))
  this%field_names = (/"VTEC"/)

  this%tme_data_file = construct_reg_grid_dataset_root(nc_file=config%tme_grid_file)
  ! no altmin/altmax: VTEC has no vertical structure, so the file's "alt"
  ! dimension is expected to have a single (placeholder) entry, and the
  ! default alt_first=1, alt_last=size(alt) covers exactly that one level.

  this%tme_obs = construct_reg_grid_dataset_group(this%tme_data_file, &
                                        field_names=this%field_names, &
                                        grid_name=trim(config%name)//" observation")
  this%tme_std = construct_reg_grid_dataset_group(this%tme_data_file, &
                                        field_names=(/"VTEC_STD"/), &
                                        grid_name=trim(config%name)//" standard deviation")

  call this%spatial_domain%init(domain_name=trim(config%name),&
                                save_members=cfg_output%save_members,&
                                kmax=cfg_output%max_moment,&
                                write_every_sec=config%write_every_sec,&
                                force_write_on_update=config%force_write_on_update,&
                                save_n_steps_after_update=cfg_output%save_n_steps_after_update)

  call this%spatial_domain%link_dataset(this%tme_data_file)
  call this%dst%init(this%tme_data_file,this%spatial_domain,this%name)
  call this%link_to_writer(this%spatial_domain, this%field_names)

  if (rank_world == 0) then
    istat = nf90_create(path='results_'//trim(config%name)//'_observations.nc',&
                      cmode=NF90_NETCDF4,&
                      ncid=this%obs_out_ncid)
    call add_global_meta_data(this%obs_out_ncid)
  end if

  if ((this%assimilate) .or. (config%always_save)) then
    ! the interpolation is already done in collect state, so that all processes can do it instead of
    ! sequential execution here on filter prcocesses
    allocate( this%computed_observations(this%tme_obs%map%size_R(mytid), n_modeltasks ) )
  end if

end subroutine obs_tum_vtec_init

!> Destroys this TUM electron-density observation's datasets and output file, and frees its arrays.
subroutine obs_tum_vtec_deallocate(this)

  ! extern
  use netcdf

  ! intern
  use mod_parallel_pdaf, only: rank_world

  implicit none

  class(obs_tum_vtec) :: this

  integer :: istat

  call this%tme_obs%destroy
  call this%tme_std%destroy
  call this%tme_data_file%destroy
  if (allocated (this%computed_observations))   deallocate (this%computed_observations)
  if (allocated (this%field_names))   deallocate (this%field_names)
  if (rank_world == 0) then
    istat = nf90_close(this%obs_out_ncid)
  end if

  call this%dst%deallocate()

end subroutine obs_tum_vtec_deallocate

!> Implementation of the observation operator for VTEC.
!!
!! Same as obs_tme_grid_obs_op_computation, but interpolates horizontally
!! only, since VTEC has no vertical extent. Consequently no geometric height
!! is required and the zg_* arguments are unused. They are only present to
!! match the interface of obs_tme_grid.
subroutine obs_tum_vtec_obs_op_computation(this, state_p, state_map, regridded, zg_mid, zg_mid_nm, zg_int, zg_int_nm)

  ! extern
  use ESMF

  ! tie-gcm
  use fields_module,only: levd0,levd1,lond0,lond1,latd0,latd1

  ! intern
  use configuration, only: cfg_filter
  use grid_observation_module, only: field_bundle, construct_field_bundle
  use quantity_info_module, only: quantity_info, get_info
  use state_module, only:  state_vector_mapping
  use tiegcm_optimized_interpolator, only: reg_grid_spline_horizontal_interpolator, &
                                           fill_halo_and_rim_2d, lon_p_halo, lat_p_halo

  implicit none

  ! ATTENTION remember to deconstruct regridded after calling this function

  class(obs_tum_vtec) :: this
  real, dimension(:), intent(in) :: state_p
  type(state_vector_mapping), intent(in) :: state_map
  type( field_bundle ), intent(inout) :: regridded
  real, dimension(levd0:levd1,lond0:lond1,latd0:latd1), intent(inout), optional :: zg_mid, zg_mid_nm, zg_int, zg_int_nm

  ! local
  type(reg_grid_spline_horizontal_interpolator) :: interp
  integer :: i
  integer :: rc
  type(quantity_info), pointer :: quantity
  type(ESMF_Field),pointer :: field

  real(ESMF_KIND_R8), dimension(:), pointer :: lon_dst, lat_dst
  real(ESMF_KIND_R8), dimension(:,:,:), contiguous, pointer :: dst_ptr ! (lon,lat,1)

  regridded = construct_field_bundle(this%tme_obs%root%grid, this%field_names)

  call ESMF_GridGetCoord(this%tme_obs%root%grid, coordDim=1, farrayPtr=lon_dst, rc=rc)
  call ESMF_GridGetCoord(this%tme_obs%root%grid, coordDim=2, farrayPtr=lat_dst, rc=rc)

  call interp%init(lon_p_halo, lat_p_halo, lon_dst, lat_dst, cfg_filter%spline_degree)

  ! loop over all fields
  do i = lbound(this%field_names,dim=1), ubound(this%field_names,dim=1)
    quantity=>get_info(this%field_names(i))
    call quantity%calc(state_map, state_p)

    ! quantities without vertical extent (LEVEL_NONE) are stored with a single
    ! level, which is exactly the shape expected by the interpolator
    call fill_halo_and_rim_2d(quantity%data(levd0,:,:))

    field=>regridded%field(name=this%field_names(i))
    call ESMF_FieldGet(field=field, localDe=0, farrayPtr=dst_ptr, rc=rc)

    call interp%interpolate(quantity%data,dst_ptr)
  end do

  call interp%destroy()

end subroutine obs_tum_vtec_obs_op_computation

!> Initialize information on the module-type observation
!!
!! The routine is called by each filter process.
!! at the beginning of the analysis step before
!! the loop through all local analysis domains.
!!
!! It has to count the number of observations of the
!! observation type handled in this module according
!! to the current time step for all observations
!! required for the analyses in the loop over all local
!! analysis domains on the PE-local state domain.
!!
!! The following four variables have to be initialized in this routine
!! * full_obs\%doassim     - Whether to assimilate this type of observations
!! * full_obs\%disttype    - type of distance computation for localization with this observaton
!! * full_obs\%ncoord      - number of coordinates used for distance computation
!! * full_obs\%id_obs_p    - index of module-type observation in PE-local state vector
!!
!! Optional is the use of
!! * full_obs\%icoeff_p    - Interpolation coefficients for obs. operator (only if interpolation is used)
!! * full_obs\%domainsize  - Size of domain for periodicity for disttype=1 (<0 for no periodicity)
!! * full_obs\%obs_err_type - Type of observation errors for particle filter and NETF (default: 0=Gaussian)
!! * full_obs\%use_global obs - Whether to use global observations or restrict the observations to the relevant ones
!!                          (default: 1=use global full observations)
!!
!! Further variables are set when the routine PDAFomi_gather_obs is called.
!!
!! **Adapting the template**
!! In this routine the variables listed above have to be initialized. One
!! can include modules from the model with 'use', e.g. for mesh information.
!! Alternatively one could include these as subroutine arguments
!!
  SUBROUTINE obs_tum_vtec_init_dim_obs(this, step, dim_obs)

    ! extern
    use esmf

    ! tie-gcm
    use mpi_module, only: mytid
    use params_module, only : nlon

    ! pdaf
    use PDAF, only: PDAFomi_gather_obs

    ! intern
    use array_print_module, only: printMat
    use array_mapping_module, only: flatten, map_to_3d
!     use mod_assimilation, only: filtertype
    use mod_assimilation, only: cutoff_radius, analysis_step_count
    use result_file_writer_module, only: instant_reg_grid_write

    IMPLICIT NONE

! *** Arguments ***
    class(obs_tum_vtec) :: this
    INTEGER, INTENT(in)    :: step       !< Current time step
    INTEGER, INTENT(inout) :: dim_obs    !< Dimension of full observation vector

! *** Local variables ***
    INTEGER :: i
    INTEGER :: dim_obs_p                 ! Number of process-local observations
    REAL, ALLOCATABLE :: obs_p(:)        ! PE-local observation vector
    REAL, ALLOCATABLE :: ivar_obs_p(:)   ! PE-local inverse observation error variance
    REAL, ALLOCATABLE :: ocoord_p(:,:)   ! PE-local observation coordinates

    real(ESMF_KIND_R8), pointer :: field_data(:,:,:)
    integer :: rc

    real, dimension(:,:,:), contiguous, pointer :: mat3d

!     real, allocatable ::debug(:,:)


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

  call this%init_full_observation_dimension

! **********************************
! *** Read PE-local observations ***
! **********************************

  ! read observation values and their coordinates
  ! also read observation error information if available
  call this%tme_obs%read_current_epoch
  call this%tme_std%read_current_epoch

! ***********************************************************
! *** Count available observations for the process domain ***
! *** and initialize index and coordinate arrays.         ***
! ***********************************************************

  ! *** Count valid observations that lie within the process sub-domain ***

  dim_obs_p = this%tme_obs%map%size_R(mytid)

  ! *** Initialize vector of observations on the process sub-domain ***

  ALLOCATE(obs_p(dim_obs_p))

  call ESMF_FieldGet(this%tme_obs%bundle%field(name='VTEC'), localDe=0, farrayPtr=field_data, &
          rc=rc)
  i=1
  obs_p( this%tme_obs%map%idx_R(i,mytid)%begin_p : &
         this%tme_obs%map%idx_R(i,mytid)%back_p ) &
        = RESHAPE( field_data, (/this%tme_obs%map%F_R_size(i,mytid)/) )

  ! VTEC is already given in tecu, matching the VTEC quantity's unit -- no conversion needed

  ! *** Initialize coordinate array of observations on the process sub-domain ***
  call this%init_coordinate_arrays(dim_obs_p, ocoord_p)

  ! *** Initialize process local index array                         ***
  ! *** This array holds the information which elements of the state ***
  ! *** vector are used in the observation operator.                 ***
  ! observations are interpolated (horizontally only, since VTEC has no
  ! vertical structure), so no direct grid-point index array is required here.
  ALLOCATE(this%full_obs%id_obs_p(0,0))

! ****************************************************************
! *** Define observation errors for process-local observations ***
! ****************************************************************

  ALLOCATE(ivar_obs_p(dim_obs_p))

  call ESMF_FieldGet(this%tme_std%bundle%field(name="VTEC_STD"), localDe=0, farrayPtr=field_data, &
          rc=rc)

  ivar_obs_p(:) = RESHAPE( field_data, (/size(field_data)/) );

  ! apply weight
  ivar_obs_p = ivar_obs_p/this%config%weight

  where (obs_p< 0.1 )
    ! samllest TIE-GCM electron density is 3100 cm-3. Thus TEC cannot be below
    ! Assuming a vertical extend of 400 km of the TIE-GCM grid
    ! 3100 cm-3 * 40 000 000 cm 1E-12 tecu/cm-2 = 0.124 tecu
    obs_p = 0.1
    ivar_obs_p = 1E+7
  end where

  ! write std
  call map_to_3d(flat=ivar_obs_p,&
                  mat3d=mat3d,&
                  mat_ub=(/this%tme_obs%root%nlon,this%tme_obs%root%nlat,this%tme_obs%root%nalt/),&
                  first=this%tme_obs%map%idx_R(i,mytid)%begin_p,&
                  last=this%tme_obs%map%idx_R(i,mytid)%back_p)

  call instant_reg_grid_write(this%obs_out_ncid,&
                              'std_'//this%tme_obs%bundle%names(i),&
                              analysis_step_count,&
                              mat3d,&
                              this%tme_data_file)
  ! write obs
  call map_to_3d(flat=obs_p,&
                  mat3d=mat3d,&
                  mat_ub=(/this%tme_obs%root%nlon,this%tme_obs%root%nlat,this%tme_obs%root%nalt/),&
                  first=this%tme_obs%map%idx_R(i,mytid)%begin_p,&
                  last=this%tme_obs%map%idx_R(i,mytid)%back_p)

  call instant_reg_grid_write(this%obs_out_ncid,&
                              this%tme_obs%bundle%names(i),&
                              analysis_step_count,&
                              mat3d,&
                              this%tme_data_file)

  ! std to inverse variance
  ivar_obs_p = 1.0/(ivar_obs_p**2);

! ****************************************
! *** Gather global observation arrays ***
! ****************************************

  CALL PDAFomi_gather_obs(this%full_obs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
       this%full_obs%ncoord, maxval(cutoff_radius(1:this%full_obs%ncoord)), dim_obs)

! ********************
! *** Finishing up ***
! ********************

  ! Deallocate all local arrays
  DEALLOCATE(obs_p, ocoord_p, ivar_obs_p)

  ! Arrays in THISOBS have to be deallocated after the analysis step
  ! by a call to deallocate_obs() in prepoststep_pdaf.

  END SUBROUTINE obs_tum_vtec_init_dim_obs

END MODULE obs_tum_vtec_pdafomi
