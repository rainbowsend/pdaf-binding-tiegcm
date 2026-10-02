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
! Defines the model_parameter type and the ~30 TIE-GCM parameters (F10.7, Kp-like drivers, chemistry rates, etc.) with their perturb/calibrate/overwrite handling.

module model_parameter_handling_module

! to enable calibration of a model parameter, the
! corresponding parameter must be included in
! `state_vector_mapping_fill_state_p_calibration` and
! `state_vector_mapping_distribute_calibration` routines

! intern
use bspline_module, only: bspline
use interpolation_module, only: nn_interpolator, spline_interpolator
use uset_module, only: str_len_char_uset

implicit none

  integer, parameter :: PARAM_HANDLING_NONE = 0       ! Do nothing with the parameter
  integer, parameter :: PARAM_HANDLING_PERTURB = 1    ! Perturb the parameter, requires file containing perturbations
  integer, parameter :: PARAM_HANDLING_CALIBRATE = 2  ! calibrate the parameter, requires file containing initial ensemble members
  integer, parameter :: PARAM_HANDLING_OVERWRITE = 3  ! overwrite the parameter, requires file containing values for each ensemble member
  integer, parameter :: PARAM_HANDLING_OVERWRITE_MEAN = 4 ! overwrite the parameter with mean value for all members

  character(len=9), dimension(0:4), parameter :: PARAM_HANDLING_CHAR &
    = (/"none     ","perturb  ","calibrate","overwrite","mean     "/)

  integer, parameter :: model_parameter_name_char_len = str_len_char_uset

  type :: model_parameter
    character(len=model_parameter_name_char_len) :: name = ""
    character(len=64) :: long_name = ""
    character(len=256) :: ensemble_file = ""
    integer :: handling = PARAM_HANDLING_NONE
    real, allocatable, dimension(:) :: val
    real, allocatable, dimension(:) :: lower_boundary
    real, allocatable, dimension(:) :: upper_boundary
    character(len=8) :: units = ""
    type(spline_interpolator), allocatable, dimension(:) :: interp
    integer :: size
    character(len=8) :: nc_dimension = "" ! only required if size > 1
    logical, dimension(1:4) :: supports_handling
    logical :: time_variable_perturbation
    contains
    procedure, pass :: allocate => model_parameter_allocate
    procedure, pass :: deallocate => model_parameter_deallocate
    procedure, pass :: update => model_parameter_update
    procedure, pass :: handle => model_parameter_handle
    procedure, pass :: constrain => model_parameter_constrain
  end type

  integer, parameter :: n_parameters = 30
  type(model_parameter), dimension(n_parameters), target :: model_parameters

  type(model_parameter), pointer, protected :: &
    param_f107, &
    param_ctpoten, &
    param_hspower, &
    param_tlbc, &
    param_zlbc, &
    param_ulbc, &
    param_vlbc, &
    param_gswm_delay, &
    param_alfac, &
    param_alfad, &
    param_colfac, &
    param_joulefac, &
    param_swden, &
    param_swvel, &
    param_imfbx, &
    param_imfby, &
    param_imfbz, &
    param_imfbf, &
!     param_igrfbu, &
!     param_igrfbe, &
!     param_igrfbn, &
!     param_igrfbf, &
!     param_igrfbd, &
!     param_igrfbi, &
    param_igrf_sh, &
    param_euvafac, &
    param_beta1, &
    param_beta2, &
    param_beta3, &
    param_beta4, &
    param_beta5, &
    param_beta6, &
    param_beta7, &
    param_beta8, &
    param_beta9, &
    param_co2u

  character(len=model_parameter_name_char_len) , dimension(:), allocatable, protected :: perturbated_parameter_names
  character(len=model_parameter_name_char_len) , dimension(:), allocatable, protected :: calibrated_parameter_names

  integer, protected :: size_dynamics = 0

contains

  !> Converts a parameter-handling name ("none"/"perturb"/"calibrate"/"overwrite"/"mean") to its PARAM_HANDLING_* integer code.
  function handling_char_to_int(handling_char) result(handling_int)
    implicit none
    character(len=9), intent(in) :: handling_char

    integer :: handling_int

    select case(handling_char)
      case('none')
        handling_int = PARAM_HANDLING_NONE
      case('perturb')
        handling_int = PARAM_HANDLING_PERTURB
      case('calibrate')
        handling_int = PARAM_HANDLING_CALIBRATE
      case('overwrite')
        handling_int = PARAM_HANDLING_OVERWRITE
      case('mean')
        handling_int = PARAM_HANDLING_OVERWRITE_MEAN
      case default
        handling_int = PARAM_HANDLING_NONE
        write(*,*) "ERROR ", handling_char, " is not a support parameter handling type"
    end select
  end function

!   subroutine get_f107(this,v)
!     use input_module,only: f107
!
!     implicit none
!     class(model_parameter) :: this
!     real, dimension(..) :: v
!
!     select rank(v)
!       rank(0)
!       v = f107
!     end select
!
!   end subroutine

  !> Defines all n_parameters TIE-GCM driver/chemistry parameters (name, units, bounds, supported handling modes) from configuration, then validates and derives the calibration parameter set.
  subroutine init_model_parameter_handling_module()

    ! intern
    use configuration, only: cfg_parameters, cfg_calibration

    implicit none
    integer :: i

    logical, dimension(n_parameters) :: calibration_mask
    logical, dimension(n_parameters) :: perturbation_mask

    character(len=128) :: msg

    write(*,*) "initalize model parameter module. Number of total parameters:", n_parameters

    i=1

    ! gpi.F90
    model_parameters(i)%name = 'f107'
    model_parameters(i)%long_name = 'F10.7 index'
    model_parameters(i)%ensemble_file = cfg_parameters%f107%ensemble_file
    call model_parameters(i)%allocate(1)
    model_parameters(i)%lower_boundary = 20.
    model_parameters(i)%upper_boundary = 600.
    model_parameters(i)%units = "sfu"
    model_parameters(i)%handling = handling_char_to_int(cfg_parameters%f107%handling)
!     model_parameters(i)%getter => get_f107
    model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .false.
    model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .true.
    param_f107 => model_parameters(i)
    i=i+1

    ! gpi.F
    model_parameters(i)%name = 'ctpoten'
    model_parameters(i)%long_name = 'cross tail potential'
    model_parameters(i)%ensemble_file = cfg_parameters%ctpoten%ensemble_file
    call model_parameters(i)%allocate(1)
    model_parameters(i)%lower_boundary = 1.
    model_parameters(i)%units = "V"
    model_parameters(i)%handling = handling_char_to_int(cfg_parameters%ctpoten%handling)
    model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .false.
    model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .true.
    param_ctpoten => model_parameters(i)
    i=i+1

    ! gpi.F
    model_parameters(i)%name = 'hspower'
    model_parameters(i)%long_name = 'hemispheric power'
    model_parameters(i)%ensemble_file = cfg_parameters%hspower%ensemble_file
    call model_parameters(i)%allocate(1)
    model_parameters(i)%lower_boundary = 1.
    model_parameters(i)%units = "GW"
    model_parameters(i)%handling = handling_char_to_int(cfg_parameters%hspower%handling)
    model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .false.
    model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .true.
    param_hspower => model_parameters(i)
    i=i+1

    ! lbc.F
    model_parameters(i)%name = 'tlbc'
    model_parameters(i)%long_name = 'lower boundary temperature'
    model_parameters(i)%ensemble_file = cfg_parameters%tlbc%ensemble_file
    call model_parameters(i)%allocate(1)
    model_parameters(i)%lower_boundary = 50.
    model_parameters(i)%units = "K"
    model_parameters(i)%handling = handling_char_to_int(cfg_parameters%tlbc%handling)
    model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .false.
    model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .false.
    param_tlbc => model_parameters(i)
    i=i+1

    ! lbc.F
    model_parameters(i)%name = 'zlbc'
    model_parameters(i)%long_name = 'lower boundary altitude'
    model_parameters(i)%ensemble_file = cfg_parameters%zlbc%ensemble_file
    call model_parameters(i)%allocate(1)
    model_parameters(i)%lower_boundary = 9500000.
    model_parameters(i)%units = "cm"
    model_parameters(i)%handling = handling_char_to_int(cfg_parameters%zlbc%handling)
    model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .false.
    model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .false.
    param_zlbc => model_parameters(i)
    i=i+1

    ! lbc.F
    model_parameters(i)%name = 'ulbc'
    model_parameters(i)%long_name = 'lower boundary zonal wind'
    model_parameters(i)%ensemble_file = cfg_parameters%ulbc%ensemble_file
    call model_parameters(i)%allocate(1)
    model_parameters(i)%units = "cm/s"
    model_parameters(i)%handling = handling_char_to_int(cfg_parameters%ulbc%handling)
    model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .false.
    model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .false.
    param_ulbc => model_parameters(i)
    i=i+1

    ! lbc.F
    model_parameters(i)%name = 'vlbc'
    model_parameters(i)%long_name = 'lower boundary meridional wind'
    model_parameters(i)%ensemble_file = cfg_parameters%vlbc%ensemble_file
    call model_parameters(i)%allocate(1)
    model_parameters(i)%units = "cm/s"
    model_parameters(i)%handling = handling_char_to_int(cfg_parameters%vlbc%handling)
    model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .false.
    model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .false.
    param_vlbc => model_parameters(i)
    i=i+1

    ! gswm.F
    model_parameters(i)%name = 'gswm_delay'
    model_parameters(i)%long_name = 'temporal delay of gswm tides'
    model_parameters(i)%ensemble_file = cfg_parameters%gswm_delay%ensemble_file
    call model_parameters(i)%allocate(1)
    model_parameters(i)%units = "s"
    model_parameters(i)%handling = handling_char_to_int(cfg_parameters%gswm_delay%handling)
    model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .false.
    model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .false.
    param_gswm_delay => model_parameters(i)
    i=i+1

    ! aurora.F
    model_parameters(i)%name = 'alfac'
    model_parameters(i)%long_name = 'Characteristic Maxwellian energy of polar cusp electrons'
    model_parameters(i)%ensemble_file = cfg_parameters%alfac%ensemble_file
    call model_parameters(i)%allocate(1)
    model_parameters(i)%lower_boundary = 0.01
    model_parameters(i)%units = "keV"
    model_parameters(i)%handling = handling_char_to_int(cfg_parameters%alfac%handling)
    model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .true.
    param_alfac => model_parameters(i)
    i=i+1

    ! aurora.F
    model_parameters(i)%name = 'alfad'
    model_parameters(i)%long_name = 'Characteristic Maxwellian energy of drizzle electrons'
    model_parameters(i)%ensemble_file = cfg_parameters%alfad%ensemble_file
    call model_parameters(i)%allocate(1)
    model_parameters(i)%lower_boundary = 0.001
    model_parameters(i)%units = "keV"
    model_parameters(i)%handling = handling_char_to_int(cfg_parameters%alfad%handling)
    model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .true.
    param_alfad => model_parameters(i)
    i=i+1

    ! input.F
    model_parameters(i)%name = 'colfac'
    model_parameters(i)%long_name = 'ion/neutral collision factor'
    model_parameters(i)%ensemble_file = cfg_parameters%colfac%ensemble_file
    call model_parameters(i)%allocate(1)
    model_parameters(i)%lower_boundary = 0.001
    model_parameters(i)%units = ""
    model_parameters(i)%handling = handling_char_to_int(cfg_parameters%colfac%handling)
    model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .true.
    param_colfac => model_parameters(i)
    i=i+1

    ! input.F
    model_parameters(i)%name = 'joulefac'
    model_parameters(i)%long_name = 'joule heating factor'
    model_parameters(i)%ensemble_file = cfg_parameters%joulefac%ensemble_file
    call model_parameters(i)%allocate(1)
    model_parameters(i)%lower_boundary = 0.001
    model_parameters(i)%units = ""
    model_parameters(i)%handling = handling_char_to_int(cfg_parameters%joulefac%handling)
    model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .true.
    param_joulefac => model_parameters(i)
    i=i+1

    ! imf.F
    model_parameters(i)%name = 'swden'
    model_parameters(i)%long_name = 'solar wind density'
    model_parameters(i)%ensemble_file = cfg_parameters%swden%ensemble_file
    call model_parameters(i)%allocate(1)
    model_parameters(i)%lower_boundary = 0.1
    model_parameters(i)%units = "1/cm3"
    model_parameters(i)%handling = handling_char_to_int(cfg_parameters%swden%handling)
    model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .false.
    model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .true.
    param_swden => model_parameters(i)
    i=i+1

    ! imf.F
    model_parameters(i)%name = 'swvel'
    model_parameters(i)%long_name = 'solar wind velocity'
    model_parameters(i)%ensemble_file = cfg_parameters%swvel%ensemble_file
    call model_parameters(i)%allocate(1)
    model_parameters(i)%lower_boundary = 50
    model_parameters(i)%units = "km/s"
    model_parameters(i)%handling = handling_char_to_int(cfg_parameters%swvel%handling)
    model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .false.
    model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .true.
    param_swvel => model_parameters(i)
    i=i+1

    ! imf.F
    model_parameters(i)%name = 'imfbx'
    model_parameters(i)%long_name = 'interplanetary magnetic field Bx at 1 AU'
    model_parameters(i)%ensemble_file = cfg_parameters%imfbx%ensemble_file
    call model_parameters(i)%allocate(1)
    model_parameters(i)%units = "nT"
    model_parameters(i)%handling = handling_char_to_int(cfg_parameters%imfbx%handling)
    model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .false.
    model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .true.
    param_imfbx => model_parameters(i)
    i=i+1

    ! imf.F
    model_parameters(i)%name = 'imfby'
    model_parameters(i)%long_name = 'interplanetary magnetic field By at 1 AU'
    model_parameters(i)%ensemble_file = cfg_parameters%imfby%ensemble_file
    call model_parameters(i)%allocate(1)
    model_parameters(i)%units = "nT"
    model_parameters(i)%handling = handling_char_to_int(cfg_parameters%imfby%handling)
    model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .false.
    model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .true.
    param_imfby => model_parameters(i)
    i=i+1

    ! imf.F
    model_parameters(i)%name = 'imfbz'
    model_parameters(i)%long_name = 'interplanetary magnetic field Bz at 1 AU'
    model_parameters(i)%ensemble_file = cfg_parameters%imfbz%ensemble_file
    call model_parameters(i)%allocate(1)
    model_parameters(i)%units = "nT"
    model_parameters(i)%handling = handling_char_to_int(cfg_parameters%imfbz%handling)
    model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .false.
    model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .true.
    param_imfbz => model_parameters(i)
    i=i+1

    ! imf.F
    model_parameters(i)%name = 'imfbf'
    model_parameters(i)%long_name = 'interplanetary magnetic total field intensity at 1 AU'
    model_parameters(i)%ensemble_file = cfg_parameters%imfbf%ensemble_file
    call model_parameters(i)%allocate(1)
    model_parameters(i)%units = "nT"
    model_parameters(i)%handling = handling_char_to_int(cfg_parameters%imfbf%handling)
    model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .false.
    model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .false.
    param_imfbf => model_parameters(i)
    i=i+1

!     ! apex.F
!     model_parameters(i)%name = 'igrfbu'
!     model_parameters(i)%long_name = 'igrf up component'
!     model_parameters(i)%ensemble_file = cfg_parameters%igrfbu%ensemble_file
!     call model_parameters(i)%allocate(1)
!     model_parameters(i)%units = "Gauss"
!     model_parameters(i)%handling = handling_char_to_int(cfg_parameters%igrfbu%handling)
!     model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .false.
!     model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .false.
!     model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .false.
!     param_igrfbu => model_parameters(i)
!     i=i+1
!
!     ! apex.F
!     model_parameters(i)%name = 'igrfbe'
!     model_parameters(i)%long_name = 'igrf east component'
!     model_parameters(i)%ensemble_file = cfg_parameters%igrfbe%ensemble_file
!     call model_parameters(i)%allocate(1)
!     model_parameters(i)%units = "Gauss"
!     model_parameters(i)%handling = handling_char_to_int(cfg_parameters%igrfbe%handling)
!     model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .false.
!     model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .false.
!     model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .false.
!     param_igrfbe => model_parameters(i)
!     i=i+1
!
!     ! apex.F
!     model_parameters(i)%name = 'igrfbn'
!     model_parameters(i)%long_name = 'igrf north component'
!     model_parameters(i)%ensemble_file = cfg_parameters%igrfbn%ensemble_file
!     call model_parameters(i)%allocate(1)
!     model_parameters(i)%units = "Gauss"
!     model_parameters(i)%handling = handling_char_to_int(cfg_parameters%igrfbn%handling)
!     model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .false.
!     model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .false.
!     model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .false.
!     param_igrfbn => model_parameters(i)
!     i=i+1
!
!     ! apex.F
!     model_parameters(i)%name = 'igrfbf'
!     model_parameters(i)%long_name = 'igrf Total Field Intensity'
!     model_parameters(i)%ensemble_file = cfg_parameters%igrfbf%ensemble_file
!     call model_parameters(i)%allocate(1)
!     model_parameters(i)%lower_boundary = 1E-16
!     model_parameters(i)%upper_boundary = 3.
!     model_parameters(i)%units = "Gauss"
!     model_parameters(i)%handling = handling_char_to_int(cfg_parameters%igrfbf%handling)
!     model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .true.
!     model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .false.
!     model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .false.
!     param_igrfbf => model_parameters(i)
!     i=i+1
!
!     ! apex.F
!     model_parameters(i)%name = 'igrfbi'
!     model_parameters(i)%long_name = 'igrf inclination'
!     model_parameters(i)%ensemble_file = cfg_parameters%igrfbi%ensemble_file
!     call model_parameters(i)%allocate(1)
!     model_parameters(i)%units = "degrees"
!     model_parameters(i)%handling = handling_char_to_int(cfg_parameters%igrfbi%handling)
!     model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .true.
!     model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .false.
!     model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .false.
!     param_igrfbi => model_parameters(i)
!     i=i+1
!
!     ! apex.F
!     model_parameters(i)%name = 'igrfbd'
!     model_parameters(i)%long_name = 'igrf declination'
!     model_parameters(i)%ensemble_file = cfg_parameters%igrfbd%ensemble_file
!     call model_parameters(i)%allocate(1)
!     model_parameters(i)%units = "degrees"
!     model_parameters(i)%handling = handling_char_to_int(cfg_parameters%igrfbd%handling)
!     model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .true.
!     model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .false.
!     model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .false.
!     param_igrfbd => model_parameters(i)
!     i=i+1

    ! apex.F
    model_parameters(i)%name = 'igrf_sh'
    model_parameters(i)%long_name = 'igrf spherical harmonic coefficients'
    model_parameters(i)%ensemble_file = cfg_parameters%igrf_sh%ensemble_file
    call model_parameters(i)%allocate(195)
    model_parameters(i)%units = "nT"
    model_parameters(i)%handling = handling_char_to_int(cfg_parameters%igrf_sh%handling)
    model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .false.
    model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .false.
    param_igrf_sh => model_parameters(i)
    i=i+1

    !qrj.F
    model_parameters(i)%name = 'euvafac'
    model_parameters(i)%long_name = 'A factor of EUVAC model'
    model_parameters(i)%ensemble_file = cfg_parameters%euvafac%ensemble_file
    call model_parameters(i)%allocate(37)
    model_parameters(i)%lower_boundary = 0.0001
    model_parameters(i)%upper_boundary = 1.
    model_parameters(i)%nc_dimension = 'bin'
    model_parameters(i)%handling = handling_char_to_int(cfg_parameters%euvafac%handling)
    model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .false.
    model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .false.
    param_euvafac => model_parameters(i)
    i=i+1

    ! chemrates_module
    model_parameters(i)%name = 'beta1'
    model_parameters(i)%long_name = 'N4S + O2 -> NO  + O + 1.4  eV'
    model_parameters(i)%ensemble_file = cfg_parameters%beta1%ensemble_file
    call model_parameters(i)%allocate(1)
    model_parameters(i)%lower_boundary = 1.50e-12
    model_parameters(i)%upper_boundary = 1.50e-10
    model_parameters(i)%handling = handling_char_to_int(cfg_parameters%beta1%handling)
    model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .false.
    model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .false.
    param_beta1 => model_parameters(i)
    i=i+1

    ! chemrates_module
    ! TODO currently no longer supported due to change in TIEGCM3
    model_parameters(i)%name = 'beta2'
    model_parameters(i)%long_name = 'N2D + O2 -> NO  + O1D + 1.84 eV'
    model_parameters(i)%ensemble_file = cfg_parameters%beta2%ensemble_file
    call model_parameters(i)%allocate(1)
    model_parameters(i)%lower_boundary = 5.00e-13
    model_parameters(i)%upper_boundary = 5.00e-11
    model_parameters(i)%handling = handling_char_to_int(cfg_parameters%beta2%handling)
    model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .false.
    model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .false.
    model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .false.
    param_beta2 => model_parameters(i)
    i=i+1

    ! chemrates_module
    model_parameters(i)%name = 'beta3'
    model_parameters(i)%long_name = 'N4S + NO -> N2  + O + 3.25 eV'
    model_parameters(i)%ensemble_file = cfg_parameters%beta3%ensemble_file
    call model_parameters(i)%allocate(1)
    model_parameters(i)%lower_boundary = 3.40e-12
    model_parameters(i)%upper_boundary = 3.40e-10
    model_parameters(i)%handling = handling_char_to_int(cfg_parameters%beta3%handling)
    model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .false.
    model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .false.
    param_beta3 => model_parameters(i)
    i=i+1

    ! chemrates_module
    model_parameters(i)%name = 'beta4'
    model_parameters(i)%long_name = 'N2D + O  -> N4S + O   + 2.38 eV'
    model_parameters(i)%ensemble_file = cfg_parameters%beta4%ensemble_file
    call model_parameters(i)%allocate(1)
    model_parameters(i)%lower_boundary = 7.00e-14
    model_parameters(i)%upper_boundary = 7.00e-12
    model_parameters(i)%handling = handling_char_to_int(cfg_parameters%beta4%handling)
    model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .false.
    model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .false.
    param_beta4 => model_parameters(i)
    i=i+1

    ! chemrates_module
    model_parameters(i)%name = 'beta5'
    model_parameters(i)%long_name = 'N2D + e  -> N4S + e + 2.38 eV'
    model_parameters(i)%ensemble_file = cfg_parameters%beta5%ensemble_file
    call model_parameters(i)%allocate(1)
    model_parameters(i)%lower_boundary = 3.60e-11
    model_parameters(i)%upper_boundary = 3.60e-09
    model_parameters(i)%handling = handling_char_to_int(cfg_parameters%beta5%handling)
    model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .false.
    model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .false.
    param_beta5 => model_parameters(i)
    i=i+1

    ! chemrates_module
    model_parameters(i)%name = 'beta6'
    model_parameters(i)%long_name = 'N2D + NO -> N2  + O   + 5.63 eV'
    model_parameters(i)%ensemble_file = cfg_parameters%beta6%ensemble_file
    call model_parameters(i)%allocate(1)
    model_parameters(i)%lower_boundary = 7.00e-12
    model_parameters(i)%upper_boundary = 7.00e-10
    model_parameters(i)%handling = handling_char_to_int(cfg_parameters%beta6%handling)
    model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .false.
    model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .false.
    param_beta6 => model_parameters(i)
    i=i+1

    ! chemrates_module
    model_parameters(i)%name = 'beta7'
    model_parameters(i)%long_name = 'N2D      -> N4S + hv'
    model_parameters(i)%ensemble_file = cfg_parameters%beta7%ensemble_file
    call model_parameters(i)%allocate(1)
    model_parameters(i)%lower_boundary = 1.06e-06
    model_parameters(i)%upper_boundary = 1.06e-04
    model_parameters(i)%handling = handling_char_to_int(cfg_parameters%beta7%handling)
    model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .false.
    model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .false.
    param_beta7 => model_parameters(i)
    i=i+1

    ! chemrates_module
    model_parameters(i)%name = 'beta8'
    model_parameters(i)%long_name = 'NO  + hv -> N4S + O'
    model_parameters(i)%ensemble_file = cfg_parameters%beta8%ensemble_file
    call model_parameters(i)%allocate(1)
    model_parameters(i)%lower_boundary = 4.50e-07
    model_parameters(i)%upper_boundary = 4.50e-05
    model_parameters(i)%handling = handling_char_to_int(cfg_parameters%beta8%handling)
    model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .false.
    model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .false.
    param_beta8 => model_parameters(i)
    i=i+1

    ! chemrates_module
    model_parameters(i)%name = 'beta9'
    model_parameters(i)%long_name = 'NO  + hv/Ly-a -> NO+ + e'
    model_parameters(i)%ensemble_file = cfg_parameters%beta9%ensemble_file
    call model_parameters(i)%allocate(1)
    model_parameters(i)%lower_boundary = 2.91e+10
    model_parameters(i)%upper_boundary = 2.91e+12
    model_parameters(i)%handling = handling_char_to_int(cfg_parameters%beta9%handling)
    model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .false.
    model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .false.
    param_beta9 => model_parameters(i)
    i=i+1

    ! newton.F
    model_parameters(i)%name = 'co2u'
    model_parameters(i)%long_name = 'CO2 concentration'
    model_parameters(i)%ensemble_file = cfg_parameters%co2u%ensemble_file
    call model_parameters(i)%allocate(1)
    model_parameters(i)%lower_boundary = 200.0 ! 280 ppm before idustrial emission
    model_parameters(i)%upper_boundary = 1000.0
    model_parameters(i)%units = "ppm"
    model_parameters(i)%handling = handling_char_to_int(cfg_parameters%co2u%handling)
    model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .true.
    model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .true.
    param_co2u => model_parameters(i)
!     i=i+1

!     ! cons_module
!     model_parameters(i)%name = 'difk'
!     model_parameters(i)%long_name = 'eddy diffusion'
!     call model_parameters(i)%allocate(1)
!     model_parameters(i)%lower_boundary = 0
!     model_parameters(i)%handling = handling_char_to_int(cfg_parameters%difk%handling)
!     model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .false.
!     model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .true.
!     model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .false.
!     param_difk => model_parameters(i)
!     i=i+1
!
!     model_parameters(i)%name = 'exp_diff_fac_o2'
!     model_parameters(i)%long_name = 'diffusion exponential factor O2'
!     call model_parameters(i)%allocate(1)
!     model_parameters(i)%lower_boundary = 1.6
!     model_parameters(i)%upper_boundary = 1.75
!     model_parameters(i)%handling = handling_char_to_int(cfg_parameters%exp_diff_fac_o2%handling)
!     model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .false.
!     model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .true.
!     model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .false.
!     param_exp_diff_fac_o2 => model_parameters(i)
!     i=i+1
!
!     model_parameters(i)%name = 'exp_diff_fac_o1'
!     model_parameters(i)%long_name = 'diffusion exponential factor O1'
!     call model_parameters(i)%allocate(1)
!     model_parameters(i)%lower_boundary = 1.6
!     model_parameters(i)%upper_boundary = 1.75
!     model_parameters(i)%handling = handling_char_to_int(cfg_parameters%exp_diff_fac_o1%handling)
!     model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .false.
!     model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .true.
!     model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .false.
!     param_exp_diff_fac_o1 => model_parameters(i)
!     i=i+1
!
!     model_parameters(i)%name = 'exp_diff_fac_n2'
!     model_parameters(i)%long_name = 'diffusion exponential factor N2'
!     call model_parameters(i)%allocate(1)
!     model_parameters(i)%lower_boundary = 1.6
!     model_parameters(i)%upper_boundary = 1.75
!     model_parameters(i)%handling = handling_char_to_int(cfg_parameters%exp_diff_fac_n2%handling)
!     model_parameters(i)%supports_handling(PARAM_HANDLING_PERTURB) = .false.
!     model_parameters(i)%supports_handling(PARAM_HANDLING_CALIBRATE) = .true.
!     model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE) = .false.
!     param_exp_diff_fac_n2 => model_parameters(i)

    if(i/=n_parameters)then
      write(msg,*) 'initalized ', i , ' parameters but expected ', n_parameters
      call shutdown("model_parameters was not allocated correctly! "//msg)
    end if

    do i = 1, n_parameters
       model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE_MEAN) =&
       model_parameters(i)%supports_handling(PARAM_HANDLING_OVERWRITE)
    end do

    ! check handling is valid
    do i = 1, n_parameters
      if(model_parameters(i)%handling/=PARAM_HANDLING_NONE) then
        if( model_parameters(i)%supports_handling(model_parameters(i)%handling).eqv..false.) then
          write(*,*) "ERROR: handling '", trim(PARAM_HANDLING_CHAR(model_parameters(i)%handling)), &
                     "' is currently not supported for '", trim(model_parameters(i)%name), &
                     "' setting handling to 'none'"
          model_parameters(i)%handling = PARAM_HANDLING_NONE
        end if
      end if
    end do

    calibration_mask = model_parameters%handling == PARAM_HANDLING_CALIBRATE
    perturbation_mask = model_parameters%handling == PARAM_HANDLING_PERTURB

    if(cfg_calibration%apply .eqv. .false.) then
      where (calibration_mask) model_parameters%handling = PARAM_HANDLING_NONE
      calibration_mask = .false.
    end if

    allocate (calibrated_parameter_names,source=pack(model_parameters%name,calibration_mask))

    do i = 1, n_parameters
      if( model_parameters(i)%handling == PARAM_HANDLING_CALIBRATE) then
        size_dynamics = size_dynamics + model_parameters(i)%size
      end if
    end do

  end subroutine

  !> Deallocates every entry in model_parameters.
  subroutine deallocate_model_parameter_handling_module
    implicit none

    integer :: i

    do i=1,n_parameters
      call model_parameters(i)%deallocate()
    end do

  end subroutine

  !> Prints each parameter's name, handling mode, and (if handled) current value to stdout; only_active restricts this to parameters with handling /= none.
  subroutine write_parmeter_handling(only_active)
    implicit none

    logical, intent(in), optional :: only_active

    ! local
    integer :: i

    do i = 1, n_parameters
      if(present(only_active)) then
        if((only_active.eqv..true.) .and. (model_parameters(i)%handling == PARAM_HANDLING_NONE)) cycle
      end if
      write(*,'(5a)') '*  ', trim(model_parameters(i)%long_name), ' (', trim(model_parameters(i)%name), ')'
      write(*,'(3x,2a)') 'handling: ', PARAM_HANDLING_CHAR(model_parameters(i)%handling)
      if(model_parameters(i)%handling /= PARAM_HANDLING_NONE) then
        write(*,'(3x,a,*(g11.4),1x,a)') 'current val: ', model_parameters(i)%val, model_parameters(i)%units
      end if
      write(*,*) ''
    end do
  end subroutine

  !> (Re)allocates this parameter's value/bounds/interpolator arrays to the given size, resetting the bounds to "unset" (spval).
  subroutine model_parameter_allocate(this,size)

    ! tie-gcm
    use params_module, only: spval ! TODO replace with infinity iee_arithmetic

    implicit none

    ! arguments
    class(model_parameter) :: this
    integer, intent(in) :: size

    call this%deallocate()

    this%size = size

    allocate(this%val(this%size))
    allocate(this%lower_boundary(this%size))
    allocate(this%upper_boundary(this%size))
    allocate(this%interp(this%size))

    this%lower_boundary = spval
    this%upper_boundary = spval
  end subroutine

  !> Deallocates this parameter's value/bounds/interpolator arrays.
  subroutine model_parameter_deallocate(this)
    implicit none

    ! arguments
    class(model_parameter) :: this

    if(allocated(this%val)) deallocate(this%val)
    if(allocated(this%lower_boundary)) deallocate(this%lower_boundary)
    if(allocated(this%upper_boundary)) deallocate(this%upper_boundary)
    if(allocated(this%interp)) deallocate(this%interp)

  end subroutine

  !> Clamps a scalar, vector, or field-shaped value (rank 0/1/2) in place to this parameter's lower/upper boundaries, where set.
  subroutine model_parameter_constrain(this, value)

    ! tie-gcm
    use params_module, only: spval ! TODO replace with infinity iee_arithmetic

    implicit none
    class(model_parameter), intent(in) :: this
    real, dimension(..), intent(inout) :: value

   ! ATTENTION parameters can be scalars and vectors. This means
   ! one can add scalar perturbations to fields of arbitrary shape
   ! and vectors to vectors

    select rank (value)
      rank(0)
        if(this%lower_boundary(1) /= spval) then
          if(value<this%lower_boundary(1)) value = this%lower_boundary(1)
        end if
        if(this%upper_boundary(1) /= spval) then
          if(value>this%upper_boundary(1)) value = this%upper_boundary(1)
        end if
      rank(1)
        where( (this%lower_boundary /= spval) .and. (value<this%lower_boundary)) value = this%lower_boundary
        where( (this%upper_boundary /= spval) .and. (value>this%upper_boundary)) value = this%upper_boundary
      rank(2)
        if(this%lower_boundary(1) /= spval) then
          where( value<this%lower_boundary(1)) value = this%lower_boundary(1)
        end if
        if(this%upper_boundary(1) /= spval) then
          where( value>this%upper_boundary(1) ) value = this%upper_boundary(1)
        end if
    end select

  end subroutine model_parameter_constrain

  !> Applies this parameter's current handling.
  subroutine model_parameter_handle(this,value)
    implicit none
    class(model_parameter), intent(in) :: this
    real, dimension(..), intent(inout) :: value


   ! ATTENTION parameters can be scalars and vectors. This means
   ! one can add scalar perturbations to fields of arbitrary shape
   ! and vectors to vectors

   if((rank(value)/=1).and.(this%size>1))then
    call shutdown("parameter has incompatible size")
   end if

    select case(this%handling)
    case(PARAM_HANDLING_PERTURB)
      select rank(value)
        rank(0)
        value = value + this%val(1)
        rank(1)
        if(this%size==1)then ! scalar perturbation
          value = value + this%val(1)
        else ! this%size has to be equatl to size(value)
          value = value + this%val(1:this%size)
        end if
        rank(2)
        value = value + this%val(1)
      end select
    case(PARAM_HANDLING_OVERWRITE,PARAM_HANDLING_OVERWRITE_MEAN)
        select rank(value)
        rank(0)
        value = this%val(1)
        rank(1)
        if(this%size==1)then ! scalar perturbation
          value = this%val(1)
        else ! this%size has to be equatl to size(value)
          value = this%val(1:this%size)
        end if
        rank(2)
        value = this%val(1)
      end select
    end select
    if (this%handling/=PARAM_HANDLING_NONE) then
      call this%constrain(value)
    end if

  end subroutine

  !> If this parameter has a handling mode and time-variable perturbation, updates its value(s) by interpolating to the current model time.
  subroutine model_parameter_update(this)

    use time_module, only: get_current_modeltime

    implicit none

    ! arguments
    class(model_parameter) :: this

    ! local
    integer :: i
    real :: seconds_since_inital_time

    if(this%handling/=PARAM_HANDLING_NONE .and. this%time_variable_perturbation) then

      call get_current_modeltime(seconds_since_inital_time)

!       write(*,*) ' vvvvvvvvvvvvvvv updating parameter ', trim(this%long_name), ' vvvvvvvvvvvvvvv '
!       write(*,*) '    seconds since inital time: ', seconds_since_inital_time
      do i=1,this%size
        this%val(i) = this%interp(i)%interpolate(seconds_since_inital_time)
      end do
!       write(*,*) '    values: ', this%val
!       write(*,*) ' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ '
    end if

  end subroutine

  !> Calls update() on every non-calibration parameter (calibration parameters are only needed at initialization).
  subroutine update_model_parameters()

    implicit none

    ! local
    integer :: i

    do i=1,n_parameters
      ! calibration parameters are only required during initalization. No need to update them.
      if(model_parameters(i)%handling/=PARAM_HANDLING_CALIBRATE)then
        call model_parameters(i)%update()
      end if
    end do

  end subroutine

end module model_parameter_handling_module
