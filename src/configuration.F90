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
! Defines and reads the namelist-based assimilation configuration.

module configuration
  implicit none

  integer, parameter :: path_len = 256

  integer, parameter :: key_len = 6
  character(len=key_len), dimension(:), allocatable, protected :: assimilated_field_names

  type :: model_parameter_namelist
    character(len=9) :: handling = "none"
    character(len=path_len) :: ensemble_file = "default"
  end type

  type :: cal_den_tme_file
    logical :: apply = .false.
    logical :: always_save = .false.
    character(len=key_len), dimension(10) :: fields=""
    character(len=path_len) :: tme_grid_file = ""
    character(len=64) :: name = ""
    real :: weight = 1.0
    real :: lb_height = 90e+3
    real :: ub_height = 700e+3
    character(len=16) :: satellite = ""
    integer :: horz_weight_half_life = 0
    integer :: vert_weight_half_life = 0
    character(len=16) :: correlations = "none"
    integer :: write_every_sec = -1
    logical :: force_write_on_update = .true.
  end type

  type :: tum_ne_tme_file
    logical :: apply = .false.
    logical :: always_save = .false.
    character(len=path_len) :: tme_grid_file = ""
    character(len=64) :: name = ""
    real :: weight = 1.0
    real :: lb_height = 150e+3
    real :: ub_height = 400e+3
    real :: min_std = 10e+3
    integer :: write_every_sec = -1
    logical :: force_write_on_update = .true.
  end type

  type :: tum_vtec_file
    logical :: apply = .false.
    logical :: always_save = .false.
    character(len=path_len) :: tme_grid_file = ""
    character(len=64) :: name = ""
    real :: weight = 1.0
    integer :: write_every_sec = -1
    logical :: force_write_on_update = .true.
  end type

  type :: satellite_file
    logical :: apply = .false.
    logical :: always_save = .false.
    character(len=64) :: file_format = "igg"
    character(len=path_len) :: file = ""
    character(len=64) :: name = ""
    real :: weight = 1.0
    integer :: write_every_sec = -1
    logical :: force_write_on_update = .true.
  end type

  ! namelist file contains objects. If a variable is not specified
  ! the default value, defined in the type is used

  type :: config_output
    character(len=64) :: result_file_name_tag = "results"
    character(len=16), dimension(20) :: saved_fields = ""
    logical :: save_state = .false.
    logical :: save_state_nm = .false.
    logical :: save_obs = .true.
    logical :: save_members = .false.
    integer :: sync_every = 10
    integer :: output_strategy = 1
    integer :: max_moment = 2
    logical :: supress_tiegcm_output = .true.
    logical :: use_double_precision = .false.
    integer :: write_every_sec = -1
    logical :: force_write_on_update = .true.
    logical :: lock_nc_time = .true.
    ! For debugging
    logical :: enforce_sync_after_update = .false.
    integer :: save_n_steps_after_update = 0
    logical :: save_unconstrained_analysis = .false.
    logical :: debug = .false.
  end type

  integer, parameter :: skip_list_len = 32
  type config_parameter
    character(len=path_len) :: ensemble_file = ""
    integer, dimension(skip_list_len) :: skip_list = 0
    type(model_parameter_namelist) :: f107
    type(model_parameter_namelist) :: ctpoten
    type(model_parameter_namelist) :: hspower
    type(model_parameter_namelist) :: tlbc
    type(model_parameter_namelist) :: zlbc
    type(model_parameter_namelist) :: ulbc
    type(model_parameter_namelist) :: vlbc
    type(model_parameter_namelist) :: gswm_delay
    type(model_parameter_namelist) :: alfac
    type(model_parameter_namelist) :: alfad
    type(model_parameter_namelist) :: colfac
    type(model_parameter_namelist) :: joulefac
    type(model_parameter_namelist) :: swden
    type(model_parameter_namelist) :: swvel
    type(model_parameter_namelist) :: imfbx
    type(model_parameter_namelist) :: imfby
    type(model_parameter_namelist) :: imfbz
    type(model_parameter_namelist) :: imfbf
    type(model_parameter_namelist) :: igrfbu
    type(model_parameter_namelist) :: igrfbe
    type(model_parameter_namelist) :: igrfbn
    type(model_parameter_namelist) :: igrfbf
    type(model_parameter_namelist) :: igrfbi
    type(model_parameter_namelist) :: igrfbd
    type(model_parameter_namelist) :: igrf_sh
    type(model_parameter_namelist) :: euvafac
    type(model_parameter_namelist) :: beta1
    type(model_parameter_namelist) :: beta2
    type(model_parameter_namelist) :: beta3
    type(model_parameter_namelist) :: beta4
    type(model_parameter_namelist) :: beta5
    type(model_parameter_namelist) :: beta6
    type(model_parameter_namelist) :: beta7
    type(model_parameter_namelist) :: beta8
    type(model_parameter_namelist) :: beta9
    type(model_parameter_namelist) :: co2u
  end type

  type config_constraints
    real :: relative_gradient_limit = -1
    logical :: quasi_neutral_ionosphere = .false.
  end type

  type :: config_calibration
    logical :: apply = .false.
    real :: localization_tapering = 1.0
    integer :: every = 1
  end type

  type :: config_filter
    logical :: open_loop = .false.
    integer :: spline_degree = 3
    integer :: first_analysis_step_sec = 0
    integer :: forecast_duration_sec = 0
    integer :: filtertype = 6
    integer :: subtype = 0
    integer :: type_trans = 0
    integer :: type_forget = 0
    real :: forget = 1.0
    integer :: type_sqrt = 0
    integer :: rank_analysis_enkf = 0
    integer :: locweight = 1
    real, dimension(3) :: cutoff_radius = 0.0
    real, dimension(3) :: support_radius = 0.0
    integer :: localization_coord_sys = 1
    integer :: sub_domain_size_vertical = 3
    integer :: sub_domain_size_zonal = 3
    integer :: sub_domain_size_meridional = 3
  end type

  type :: config_state
    logical :: o2 = .false.
    logical :: o1 = .false.
    logical :: he = .false.
    logical :: tn = .false.
    logical :: ne = .false.
    logical :: zonal_wind = .false.
    logical :: meridional_wind = .false.
    logical :: atomic_oxygen_ion_density = .false.
    logical :: molecular_oxygen_ion_density = .false.
    logical :: atomic_argon = .false.
    logical :: nitric_oxide = .false.
    logical :: excited_atomic_nitrogen_4s = .false.
    logical :: excited_atomic_nitrogen_2d = .false.
    logical :: electron_temperature = .false.
    logical :: ion_temperature = .false.
    logical :: vertical_motion = .false.
    logical :: electric_potential = .false.
  end type

  type :: config_ensemble
    integer :: ensemble_size = 0
    logical :: overwrite_source = .false.
    character(len=path_len) :: source_path = ""
    character(len=path_len) :: source_name = ""
  end type

  integer, parameter :: max_cal_den_tme_file = 3
  integer, parameter :: max_satellite_file = 6

  type :: config_observation
    type(cal_den_tme_file), dimension(max_cal_den_tme_file) :: cal_den
    type(satellite_file), dimension(max_satellite_file) :: satellite
    type(tum_ne_tme_file) :: tum_ne
    type(tum_vtec_file) :: tum_vtec
  end type

  type :: config_log
      integer :: verbose_level = 0
  end type

  type(config_output), protected :: cfg_output
  type(config_calibration), protected :: cfg_calibration
  type(config_constraints), protected :: cfg_constraints
  type(config_parameter), protected :: cfg_parameters
  type(config_filter), protected :: cfg_filter
  type(config_state), protected :: cfg_state
  type(config_ensemble), protected :: cfg_ensemble
  type(config_observation), protected :: cfg_observation
  type(config_log), protected :: cfg_log

  contains

  !> Reads the namelist-based assimilation configuration file and populates
  !! the module's cfg_* variables, deriving assimilated_field_names and
  !! validating settings.
  subroutine read_config( namelistfile )

    ! tie-gcm
    use mpi_module, only: mytid

    ! args
    character(len=*), intent(in) :: namelistfile

    ! local
    integer :: dim_assimilated_fields
    logical :: THERE
    integer :: i
    logical, dimension(17) :: mask_state

    type(config_output) :: output
    type(config_calibration) :: calibration
    type(config_constraints) :: constraints
    type(config_parameter) :: parameters
    type(config_filter) :: filter
    type(config_state) :: state
    type(config_ensemble) :: ensemble
    type(config_observation) :: observation
    type(config_log) :: logger



    namelist / da_settings / &
        output, &
        calibration, &
        constraints, &
        parameters, &
        filter, &
        state, &
        ensemble, &
        observation,&
        logger


    INQUIRE( FILE=trim(namelistfile), EXIST=THERE )
    if( .not. THERE) then
        call shutdown('required assimilation configuration file does not exist: ' // trim(namelistfile) )
    end if

    open( 10, file = trim(namelistfile), status='old' )

    write(*,*) "reading configuration"
    read( 10, nml = da_settings )

    output%max_moment = min(output%max_moment,4)

    if(output%debug)then
      write(*,*) "!!! ATTENTION output%debug is enabled. This may slow down the program significantly."
      if(.not.output%save_unconstrained_analysis) then
        write(*,*) "output%debug is enabled. Enforcing output%save_unconstrained_analysis option"
        output%save_unconstrained_analysis = .true.
      end if
      if(.not.output%save_obs) then
        write(*,*) "output%debug is enabled. Enforcing output%save_obs option"
        output%save_obs = .true.
      end if
      if(.not.output%save_state) then
        write(*,*) "output%debug is enabled. Enforcing output%save_state option"
        output%save_state = .true.
      end if
      if(.not.output%save_state_nm) then
        write(*,*) "output%debug is enabled. Enforcing output%save_state_nm option"
        output%save_state_nm = .true.
      end if
      if(.not.output%save_members) then
        write(*,*) "output%debug is enabled. Enforcing output%save_members option"
        output%save_members = .true.
      end if
      if(output%save_n_steps_after_update < 5) then
        write(*,*) "output%debug is enabled. Raising output%save_n_steps_after_update option to 5"
        output%save_n_steps_after_update = 5
      end if
      if(.not.output%force_write_on_update) then
        write(*,*) "output%debug is enabled. Enforcing output%force_write_on_update option"
        output%force_write_on_update = .true.
      end if
      if(.not.output%enforce_sync_after_update) then
        write(*,*) "output%debug is enabled. Enforcing output%enforce_sync_after_update option"
        output%enforce_sync_after_update = .true.
      end if
      if( output%sync_every/=1) then
        write(*,*) "output%debug is enabled. Setting output%sync_every to 1"
        output%sync_every=1
      end if
      if(logger%verbose_level < 1) then
        write(*,*) "output%debug is enabled. Raising logger%verbose_level option to 1"
        logger%verbose_level = 1
      end if
    end if

    if(output%lock_nc_time) then
      write(*,*) "ATTENTION lock_nc_time has been activated. write_every_sec in observation settings are ignored."
      do i=1,max_cal_den_tme_file
        observation%cal_den(i)%write_every_sec= output%write_every_sec
        observation%cal_den(i)%force_write_on_update= output%force_write_on_update
      end do
      do i=1,max_satellite_file
        observation%satellite(i)%write_every_sec= output%write_every_sec
        observation%satellite(i)%force_write_on_update= output%force_write_on_update
      end do
      observation%tum_ne%write_every_sec= output%write_every_sec
      observation%tum_ne%force_write_on_update= output%force_write_on_update
    else
      do i=1,max_cal_den_tme_file
        if(observation%cal_den(i)%write_every_sec==-1 ) observation%cal_den(i)%write_every_sec= output%write_every_sec
      end do
      do i=1,max_satellite_file
        if(observation%satellite(i)%write_every_sec==-1 ) observation%satellite(i)%write_every_sec= output%write_every_sec
      end do
      if(observation%tum_ne%write_every_sec==-1 ) observation%tum_ne%write_every_sec= output%write_every_sec
    end if

    ! perform some consitency checks 
    if( filter%open_loop .eqv. .false.) then

      mask_state=(/state%o2,&
                   state%o1,&
                   state%he,&
                   state%tn,&
                   state%ne,&
                   state%zonal_wind,&
                   state%meridional_wind,&
                   state%atomic_oxygen_ion_density,&
                   state%molecular_oxygen_ion_density,&
                   state%atomic_argon,&
                   state%nitric_oxide,&
                   state%excited_atomic_nitrogen_4s,&
                   state%excited_atomic_nitrogen_2d,&
                   state%electron_temperature,&
                   state%ion_temperature,&
                   state%vertical_motion,&
                   state%electric_potential/)
      dim_assimilated_fields = count( mask_state )

      if( dim_assimilated_fields == 0) then
          write(*,*) 'WARNING! no TIE-GCM field for assimilation selected'
      end if

      allocate( assimilated_field_names(dim_assimilated_fields) )

      assimilated_field_names = pack( (/"O2    ",&
                                        "O1    ",&
                                        "HE    ",&
                                        "TN    ",&
                                        "NE    ",&
                                        "UN    ",&
                                        "VN    ",&
                                        "OP    ",&
                                        "O2P   ",&
                                        "AR    ",&
                                        "NO    ",&
                                        "N4S   ",&
                                        "N2D   ",&
                                        "TE    ",&
                                        "TI    ",&
                                        "OMEGA ",&
                                        "POTEN "/),&
                                        mask_state )

      write(*,'(a,1x,*(a4))') 'the following TIE-GCM fields are part of the state vector' , assimilated_field_names

      if( observation%tum_ne%apply .and. .not. state%ne) then
          write(*,*) 'WARNING! tum_ne is set to true, but ne is not part of state vector.'
      end if

      if( observation%tum_ne%ub_height < observation%tum_ne%lb_height ) then
        call shutdown('observation%tum_ne%ub_height must be larger than or equal to observation%tum_ne%lb_height' )
      end if

      do i=1, max_cal_den_tme_file
        if( observation%cal_den(i)%ub_height < observation%cal_den(i)%lb_height ) then
          call shutdown('observation%cal_den%ub_height must be larger than or equal to observation%cal_den%lb_height' )
        end if
      end do
    end if

    if( filter%open_loop .eqv. .false.) then
      where(filter%cutoff_radius<0) filter%cutoff_radius = maxval(filter%cutoff_radius)
      where(filter%support_radius<0) filter%support_radius = maxval(filter%support_radius)
    end if

    if(calibration%apply)then
      if(calibration%every<1) then
        call shutdown('calibration%every must be larger than zero' )
      end if
    end if

    ! assign to module objects
    cfg_output = output
    cfg_filter = filter
    cfg_ensemble = ensemble
    cfg_calibration = calibration
    cfg_parameters = parameters
    cfg_state = state
    cfg_constraints = constraints
    cfg_observation = observation
    cfg_log = logger

    if(mytid==0) write( *, nml = da_settings )

  end subroutine read_config

  !> Deallocates variables that were allocated when reading the configuration. Called at the end of the program.
  subroutine deallocate_configuration
    implicit none
    if ( allocated( assimilated_field_names) ) deallocate(assimilated_field_names)
  end subroutine deallocate_configuration

end module configuration
