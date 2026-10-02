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
! PDAF callback: sets the state dimension and coordinates of the current local analysis domain.

!BOP
!
! !ROUTINE: init_dim_l_pdaf --- Set dimension of local model state
!
! !INTERFACE:
SUBROUTINE init_dim_l_pdaf(step, domain_p, dim_l)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called during analysis step
! in the loop over all local analysis domain.
! It has to set the dimension of local model 
! state on the current analysis domain.
!
! The routine is called by each filter process.
!
! !USES:

  ! tie-gcm
  use mpi_module, only: mytid

  ! intern
  use configuration,&
      only: cfg_filter
  use cell_id_coordinate_system, only: to_zonal_meridional_vertical
  use mod_assimilation,&
      only: coords_l, analysis_step_count
  use tgcm_pdaf_omi_obs_type_module, only: COORD_SPH_3D, COORD_CELL_IDX_3D
  use quantity_computation_module, only: zg_mid_mean
!   use mod_parallel_pdaf,&
!     only: rank_filter
  use state_module,&
      only: lon_p, lat_p, idx_intern, idx_nc, levX0, levX1, state_vector
  use structured_gird_subdomain_module,&
      only: subdomains, get_center_coords, get_center_coords_cell_id

  ! pdaf
  use PDAF, only: PDAFomi_set_debug_flag, PDAFlocal_set_indices

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: step     ! Current time step
  INTEGER, INTENT(in)  :: domain_p ! Current local analysis domain
  INTEGER, INTENT(out) :: dim_l    ! Local state dimension

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_dim_l)
! Called by: PDAF_lestkf_update  (as U_init_dim_l)
! Called by: PDAF_letkf_update   (as U_init_dim_l)
! Called by: PDAF_lnetf_update   (as U_init_dim_l)
!EOP

! local
  integer, allocatable :: map(:) ! indices of local state vector elements in the PE-local global state vector
  integer :: first, last
  integer :: i

!   IF (domain_p==1600 .AND. rank_filter==0) THEN
!     CALL PDAF_set_debug_flag(domain_p)
!     CALL PDAFomi_set_debug_flag(domain_p)
!   ELSE
!     CALL PDAF_set_debug_flag(0)
!   ENDIF

! ****************************************
! *** Initialize local state dimension ***
! ****************************************

   ! Only state is localized. Global parameters are updated later in prepoststep

   dim_l = subdomains(domain_p)%domain_size * size(state_vector%tgcm_field_names)


! **********************************************
! *** Initialize coordinates of local domain ***
! **********************************************

  ! Global coordinates of local analysis domain
  ! one point represeting position of local analysis sub domain

  ! TODO could be solved nicer, without allocating each call
  if(allocated(coords_l)) DEALLOCATE(coords_l)

  select case(cfg_filter%localization_coord_sys)
      case(COORD_SPH_3D)
        ! ATTENTION zg mean is used for all members
        allocate(coords_l(3))
        ! (lon, lat, alt) tuple for current subdomain
        coords_l(:) = get_center_coords(domain_p,lon_p,lat_p,&
                                     zg_mid_mean(levX0:levX1,&
                                          idx_intern(mytid)%lon0:idx_intern(mytid)%lon1,&
                                          idx_intern(mytid)%lat0:idx_intern(mytid)%lat1))
      case(COORD_CELL_IDX_3D)
         allocate(coords_l(3))
         coords_l(:) = get_center_coords_cell_id(domain_p, idx_nc(mytid)%lon0-1, idx_nc(mytid)%lat0-1)
         ! get_center_coords_cell_id returns tuple (lev,lon,lat). to_zonal_meridional_vertical changes the
         ! tuple to (lon,lat,lev)
         call to_zonal_meridional_vertical(coords_l)
      case default
          call shutdown('init_dim_l_pdaf: invalid choice for localization coordinate system')
  end select

  ! only write once at the beginning
  if(analysis_step_count==1)then
    write(*,*) 'sub domian id:', domain_p, ' coordinates: ', coords_l
  end if

! ******************************************************
! *** Initialize array of indices of the local state ***
! ***  vector elements in the global state vector.   ***
! ******************************************************

  ! For each 3d TIE-GCM field making up the state vector, the local subdomain
  ! uses the same mapping (subdomains(domain_p)%mapping), offset by the
  ! field's own base index in the PE-local global state vector.

  allocate(map(dim_l))

  first = 1
  last = subdomains(domain_p)%domain_size
  do i = state_vector%idx_f3d_0, state_vector%idx_f3d_1
    map(first:last) = state_vector%map%idx_R(i, mytid)%begin_p - 1 + subdomains(domain_p)%mapping(:)
    first = last + 1
    last = last + subdomains(domain_p)%domain_size
  end do

  call PDAFlocal_set_indices(dim_l, map)

  deallocate(map)

END SUBROUTINE init_dim_l_pdaf
