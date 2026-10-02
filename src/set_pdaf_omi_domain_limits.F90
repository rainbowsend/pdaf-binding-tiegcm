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
! Tells PDAF-OMI the lon/lat bounding box of each MPI subdomain, needed for localized observation handling.

!> Computes and passes to PDAF-OMI the lon/lat (or cell-id-coordinate) bounding box of this MPI subdomain.
subroutine set_pdaf_omi_domain_limits()

  ! ATTENTION tiegcm_optimized_interpolator and cell_id_coordinate_system need to be initalized before calling

  ! extern
  use PDAFomi_obs_f, only: PDAFomi_set_domain_limits

  ! tiegcm
  use cons_module, only: pi
  use mpi_module, only: mytid

  ! intern
  use cell_id_coordinate_system, only: longitude_to_cell_id_coordinate, &
                                        latitude_to_cell_id_coordinate
  use configuration, only: cfg_filter
  use tgcm_pdaf_omi_obs_type_module, only: COORD_SPH_3D
  use tiegcm_optimized_interpolator, only: lb_lon, ub_lon, ub_lat, lb_lat

  implicit none

  ! local
  real, dimension(2,2) :: lim

  ! https://pdaf.awi.de/trac/wiki/OMI_use_global_obs

  select case(cfg_filter%localization_coord_sys)
    case(COORD_SPH_3D)
      lim(1,1) = lb_lon(mytid)/180.*pi ! western edge of the domain
      lim(1,2) = ub_lon(mytid)/180.*pi ! eastern edge of the domain
      lim(2,1) = ub_lat(mytid)/180.*pi ! northern edge of the domain
      lim(2,2) = lb_lat(mytid)/180.*pi ! southern edge of the domain
    case default
      lim(1,1) = longitude_to_cell_id_coordinate(lb_lon(mytid)) ! western edge of the domain
      lim(1,2) = longitude_to_cell_id_coordinate(ub_lon(mytid)) ! eastern edge of the domain
      lim(2,1) = latitude_to_cell_id_coordinate(ub_lat(mytid)) ! northern edge of the domain
      lim(2,2) = latitude_to_cell_id_coordinate(lb_lat(mytid)) ! southern edge of the domain
  end select

  write(*,*) 'set_pdaf_omi_domain_limits'
  write(*,*) 'first dim in [',lim(1,1),',',lim(1,2),']'
  write(*,*) 'second dim in [',lim(2,1),',',lim(2,2),']'

  call PDAFomi_set_domain_limits(lim)

end subroutine
