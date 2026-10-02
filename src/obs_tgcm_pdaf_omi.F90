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
! Adds TIE-GCM-specific next_observation/localization behavior shared by all TIE-GCM PDAF-OMI observation types.

! adds next_observation as deferred procedure to observation_interface
module tgcm_pdaf_omi_obs_type_module

  use pdaf_omi_obs_type_module, only: observation_interface

  implicit none

  type, extends(observation_interface), abstract :: tgcm_pdaf_observation_interface
    ! computed_observations holds the 'observations' computed from the TIE-GCM's state
    ! size: number of observations on subdomain x number of ensemble members
    real, dimension(:,:), allocatable :: computed_observations
    !
    logical :: is_available_at_next_assim_step
    ! The observation operator first computes the observation on the native TIE_GCM grid
    ! and then interpolates it.
    ! On which native TIE-GCM grid the observation is computed before interpolation
    ! to the location of the observationis saved here
    integer :: tiegcm_native_grid ! LEVEL_MID, LEVEL_INT, LEVEL_NONE
    contains
    ! derived procedures
    procedure, pass(this) :: prodRinvA_l
    ! class procedures
    procedure, pass(this) :: init_full_observation_dimension
    procedure, pass(this) :: link_to_writer
    procedure(next_observation_interface), deferred :: next_observation
    procedure(deallocate_interface), deferred :: deallocate
  end type tgcm_pdaf_observation_interface

! these procedures must be implemented by any implementation of the
! tgcm_pdaf_observation_interface
  abstract interface
    !> Deferred: reports whether this observation type has data available before next_analysis_step.
    function next_observation_interface(this, next_analysis_step) result (is_available)
      use ESMF
      import tgcm_pdaf_observation_interface
      class(tgcm_pdaf_observation_interface) :: this
      type(ESMF_Time), intent(in) :: next_analysis_step
      logical :: is_available
    end function
    !> Deferred: deallocates this observation type's own resources.
    subroutine deallocate_interface(this)
      import tgcm_pdaf_observation_interface
      class(tgcm_pdaf_observation_interface) :: this
    end subroutine
  end interface

  integer, parameter :: COORD_SPH_3D = 1
  integer, parameter :: COORD_CELL_IDX_3D = 2

  contains

  !> Default prodRinvA_l implementation (overridable in subclasses): applies the observation error covariance inverse via PDAFomi_prodRinvA_l.
  SUBROUTINE prodRinvA_l(this,domain_p, step, dim_obs_l, rank, obs_l, A_l, C_l)

    ! PDAF
    use PDAF, only: PDAFomi_prodRinvA_l
    ! the dummy argument obs_l below shadows the type name, hence the rename
    use PDAFomi_obs_l, only: obs_l_type => obs_l

    IMPLICIT NONE

    ! *** Arguments ***
    class(tgcm_pdaf_observation_interface) :: this
    INTEGER, INTENT(in) :: domain_p             !< Index of current local analysis domain
    INTEGER, INTENT(in) :: step                 !< Current time step
    INTEGER, INTENT(in) :: dim_obs_l            !< Dimension of local observation vector
    INTEGER, INTENT(in) :: rank                 !< Rank of initial covariance matrix
    REAL, INTENT(in)    :: obs_l(dim_obs_l)     !< Local vector of observations
    REAL, INTENT(inout) :: A_l(dim_obs_l, rank) !< Input matrix
    REAL, INTENT(out)   :: C_l(dim_obs_l, rank) !< Output matrix

    ! local
    type(obs_l_type), pointer :: local_obs ! this thread's local observation

    call this%local_obs(local_obs)

    CALL PDAFomi_prodRinvA_l(local_obs, this%full_obs, dim_obs_l, rank, A_l, C_l, 1)

  END SUBROUTINE prodRinvA_l

  !> Registers this observation type's fields as "observed" quantities and adds its spatial domain to the result writer.
  subroutine link_to_writer(this, spatial,field_names)

    ! intern
    use result_file_writer_module, only: spatial_domain
    use result_writer_frontend, only: result_writer
    use quantity_computation_module, only: qset

    implicit none

    ! arguments
    class(tgcm_pdaf_observation_interface) :: this
    class(spatial_domain) :: spatial
    character(len=*), dimension(:), intent(in) :: field_names

    call qset%add(field_names, "observed")
    call result_writer%add(spatial)


  end subroutine

  !> Returns the native TIE-GCM grid (LEVEL_MID, LEVEL_INT or LEVEL_NONE) the
  !! given quantity is computed on before the interpolation
  !!
  !! Deriving the value from the quantity keeps it consistent with
  !! quantity_info_module. Only usable after init_quantity_info_module,
  !! which init_pdaf calls before init_observations_module.
  function get_tiegcm_native_grid(varname) result(level)

    ! intern
    use quantity_info_module, only: quantity_info, get_info

    implicit none

    ! arguments
    character(len=*), intent(in) :: varname

    ! result
    integer :: level

    ! local
    type(quantity_info), pointer :: info

    info => get_info(varname)
    level = info%level

  end function get_tiegcm_native_grid

  !> Sets PDAF-OMI's distance-computation type, coordinate count, and domain size for this observation's localization, based on the configured coordinate system.
  subroutine init_full_observation_dimension(this)

    ! tie-gcm
    use params_module, only : nlon

    ! intern
    use configuration, only: cfg_filter
    use mod_parallel_pdaf, only: rank_filter
    use quantity_info_module, only: LEVEL_MID, LEVEL_INT, LEVEL_NONE

    implicit none

    ! arguments
    class(tgcm_pdaf_observation_interface) :: this

    IF (rank_filter==0) &
         WRITE (*,'(8x,a,a)') 'Assimilate observations - ', trim(this%name)

    ! Store whether to assimilate this observation type (used in routines below)
    IF (this%assimilate) this%full_obs%doassim = 1

    ! disttype:
    !   Specify type of distance computation
    !   see PDAF/src/PDAFomi_obs_l.F90
    !     1: periodic cartesian
    !     3: Geographic distance computation in meters using haversine formula
    !        with horizontal coordinates in radians
    !        (latitude: -pi/2 to +pi/2; longitude -pi/+pi or 0 to 2pi)
    ! ncoord:
    !   Number of coordinates used for distance computation
    !   The distance compution starts from the first row
    !
    !   Quantities without vertical extent, i.e. VTEC, are integrated over the
    !   whole column and are therefore localized horizontally only. As the
    !   vertical is the last coordinate in both systems (COORD_SPH_3D, COORD_CELL_IDX_3D),
    !   this is achieved by setting ncoord=2.
      select case(this%tiegcm_native_grid)
        case(LEVEL_MID,LEVEL_INT)
          this%full_obs%ncoord = 3
        case(LEVEL_NONE)
          this%full_obs%ncoord = 2
        case default
          call shutdown('init_full_observation_dimension: unhandled tiegcm_native_grid')
    end select


    select case(cfg_filter%localization_coord_sys)
      case(COORD_SPH_3D)
        ! order: longitude, latitude, altitude
        this%full_obs%disttype = 3
      case(COORD_CELL_IDX_3D)
        ! order: longitude, latitude, vertical. The cell indices are stored in
        ! the TIE-GCM order (vertical, longitude, latitude) and are rearranged
        ! by to_zonal_meridional_vertical, see
        ! obs_tme_grid_init_coordinate_arrays and init_dim_l_pdaf.
        this%full_obs%disttype = 1
    end select

    ! domainsize is required for distance computation (only: 1 periodic Cartesian distance)
    ! see PDAFomi_comp_dist2 in PDAFomi_obs_l.F90
    allocate(this%full_obs%domainsize(this%full_obs%ncoord))
    this%full_obs%domainsize = 0
    select case(cfg_filter%localization_coord_sys)
      case(COORD_CELL_IDX_3D)
        this%full_obs%domainsize(1) = nlon ! total number of cells in longitude (72 or 144)
    end select
  end subroutine

end module tgcm_pdaf_omi_obs_type_module
