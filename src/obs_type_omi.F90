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
! Abstract PDAF-OMI observation_interface defining the deferred procedures every observation type must implement.

module pdaf_omi_obs_type_module

  ! PDAF
  use PDAFomi_obs_f, only: obs_f
  use PDAFomi_obs_l, only: obs_l

  implicit none

  ! upper bound on simultaneously initialized observation types
  integer, parameter :: max_obs_types = 20

  ! PDAF requires the local observation (obs_l) to be private to each OpenMP
  ! thread: init_dim_obs_l and prodRinvA_l are called from inside PDAF's
  ! OpenMP loop over local analysis domains, obs_l holds per-domain scratch
  ! (id_obs_l, distance_l, ...), and since PDAF 3.1.1 PDAF deallocates its
  ! components once per domain (PDAFomi_dealloc_local). A shared obs_l is
  ! therefore raced and double-freed. See PDAF's obs_OBSTYPE_pdafomi
  ! templates, which declare THREADPRIVATE(thisobs_l).
  !
  ! obs_l cannot be a component of observation_interface: the observation
  ! objects are shared module variables, THREADPRIVATE only applies to named
  ! module/SAVE variables, and full_obs must stay shared as it is filled once
  ! outside the parallel region. Each observation type therefore owns a slot
  ! of this pool, selected by obs_l_id and accessed through the local_obs
  ! type-bound procedure below. The size must be a constant: an allocatable
  ! THREADPRIVATE array would be unallocated on all but the allocating thread.
  !
  ! The slot must be identified by an index, not by a pointer stored in the
  ! observation object: set_obs_l_id runs once on the master thread outside
  ! any parallel region, so a stored pointer would bind to the master's copy
  ! of the pool and every thread would follow it there. THREADPRIVATE
  ! duplicates the variable, not pointers aimed at it. Indexing the pool from
  ! inside local_obs is what makes the lookup resolve on the calling thread.
  type(obs_l), dimension(max_obs_types), target, private :: local_obs_pool
!$OMP THREADPRIVATE(local_obs_pool)

  ! interface definition for observations ----------------------------------------
  type, abstract :: observation_interface
    ! member variables
    type(obs_f) :: full_obs   ! full observation
    ! slot of this observation type in local_obs_pool, set once by
    ! set_obs_l_id. The obs_l itself is deliberately not a component,
    ! see local_obs_pool above.
    integer, private :: obs_l_id = 0

    logical :: assimilate   !< Whether to assimilate this data type
    real    :: rms_obs      !< Observation error standard deviation (for constant errors)

    character(len=64) :: name

    contains
    ! member procedures
    procedure, pass(this) :: set_obs_l_id
    procedure, pass(this) :: local_obs => get_local_obs
    ! arguments of deferred procedures are declared in abstract interface block
    procedure(init_dim_obs_interface), deferred :: init_dim_obs
    procedure(obs_op_interface), deferred :: obs_op
    procedure(init_dim_obs_l_interface), deferred :: init_dim_obs_l
    procedure(localize_covar_interface), deferred :: localize_covar
    procedure(prodRinvA_l_interface), deferred :: prodRinvA_l
  end type observation_interface

  ! these procedures must be implemented by any implementation of the
  ! observation_interface
  abstract interface
    subroutine init_dim_obs_interface(this,step,dim_obs)
      import observation_interface
      class(observation_interface) :: this
      INTEGER, INTENT(in)    :: step       !< Current time step
      INTEGER, INTENT(inout) :: dim_obs    !< Dimension of full observation vector
    end subroutine

    subroutine obs_op_interface(this,dim_p,dim_obs,state_p,ostate)
      import observation_interface
      class(observation_interface) :: this
      INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
      INTEGER, INTENT(in) :: dim_obs               !< Dimension of full observed state (all observed fields)
      REAL, INTENT(in)    :: state_p(dim_p)        !< PE-local model state
      REAL, INTENT(inout) :: ostate(dim_obs)       !< Full observed state
    end subroutine

    subroutine init_dim_obs_l_interface(this,domain_p,step,dim_obs,dim_obs_l)
      import observation_interface
      class(observation_interface) :: this
      INTEGER, INTENT(in)  :: domain_p     !< Index of current local analysis domain
      INTEGER, INTENT(in)  :: step         !< Current time step
      INTEGER, INTENT(in)  :: dim_obs      !< Full dimension of observation vector
      INTEGER, INTENT(inout) :: dim_obs_l  !< Local dimension of observation vector
    end subroutine

    subroutine localize_covar_interface(this,dim_p,dim_obs,HP_p,HPH,coords_p)
      import observation_interface
      class(observation_interface) :: this
      INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
      INTEGER, INTENT(in) :: dim_obs               !< Dimension of observation vector
      REAL, INTENT(inout) :: HP_p(dim_obs, dim_p)  !< PE local part of matrix HP
      REAL, INTENT(inout) :: HPH(dim_obs, dim_obs) !< Matrix HPH
      REAL, INTENT(in)    :: coords_p(:,:)         !< Coordinates of state vector elements
    end subroutine

    subroutine prodRinvA_l_interface(this,domain_p, step, dim_obs_l, rank, obs_l, A_l, C_l)
      import observation_interface
      class(observation_interface) :: this
      INTEGER, INTENT(in) :: domain_p          !< Index of current local analysis domain
      INTEGER, INTENT(in) :: step              !< Current time step
      INTEGER, INTENT(in) :: dim_obs_l         !< Dimension of local observation vector
      INTEGER, INTENT(in) :: rank              !< Rank of initial covariance matrix
      REAL, INTENT(in)    :: obs_l(dim_obs_l)  !< Local vector of observations
      REAL, INTENT(inout) :: A_l(dim_obs_l, rank) !< Input matrix
      REAL, INTENT(out)   :: C_l(dim_obs_l, rank) !< Output matrix
    end subroutine
  end interface

  contains

  !> Assigns this observation type its slot in local_obs_pool. Called once
  !! per initialized observation type by init_observations_module.
  subroutine set_obs_l_id(this, id)

    implicit none

    ! arguments
    class(observation_interface) :: this
    integer, intent(in) :: id

    if((id < 1) .or. (id > max_obs_types)) then
      call shutdown('set_obs_l_id: id outside local_obs_pool bounds')
    end if

    this%obs_l_id = id

  end subroutine

  !> Points lobs at the calling thread's copy of this observation type's
  !! local observation (obs_l).
  !!
  !! local_obs_pool is THREADPRIVATE, so indexing it here resolves to the
  !! storage of the thread executing the call. ATTENTION therefore call this
  !! inside the routine that uses the result - a pointer stored in the
  !! (shared) object or kept across parallel regions would alias one thread's
  !! copy for all threads, recreating the race this construction removes.
  !!
  !! This is deliberately a subroutine and not a pointer-valued function.
  subroutine get_local_obs(this, lobs)

    implicit none

    ! arguments
    class(observation_interface) :: this
    type(obs_l), pointer, intent(out) :: lobs

    if(this%obs_l_id < 1) then
      call shutdown('local_obs: set_obs_l_id was not called for observation '//trim(this%name))
    end if

    lobs => local_obs_pool(this%obs_l_id)

  end subroutine

end module pdaf_omi_obs_type_module
