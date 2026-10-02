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
! Computes derived diagnostic quantities (e.g. geometric height) and interpolates them onto configured observation/output grids each step.

module quantity_computation_module

! intern
use quantity_info_module, only: quantity_type
use quantity_sets_module, only: quantity_sets

use configuration, only: cfg_log

implicit none

real, allocatable, dimension(:,:,:) ::  zg_int, zg_int_nm, zg_mid, zg_mid_nm
real, allocatable, dimension(:,:,:) ::  zg_mid_mean, zg_mid_mean_nm, zg_int_mean, zg_int_mean_nm
type(quantity_sets) :: qset
type(quantity_type), dimension(:), allocatable :: quantities

contains

!> Allocates the geometric-height working arrays (zg_int/zg_mid and their
!! _nm/_mean variants) and the georeferenced output datasets.
subroutine init_quantity_computation_module()

  ! intern
  use georeferenced_data_module, only: georeferenced_data_allocate

  ! tiegcm
  use fields_module, only: levd0, levd1, lond0, lond1, latd0, latd1

  implicit none

  allocate( zg_int( levd0:levd1,lond0:lond1,latd0:latd1 ) )
  allocate( zg_mid( levd0:levd1,lond0:lond1,latd0:latd1 ) )
  allocate( zg_mid_nm( levd0:levd1,lond0:lond1,latd0:latd1 ) )
  allocate( zg_int_nm( levd0:levd1,lond0:lond1,latd0:latd1 ) )

  allocate( zg_mid_mean( levd0:levd1,lond0:lond1,latd0:latd1 ) )
  allocate( zg_mid_mean_nm( levd0:levd1,lond0:lond1,latd0:latd1 ) )
  allocate( zg_int_mean( levd0:levd1,lond0:lond1,latd0:latd1 ) )
  allocate( zg_int_mean_nm( levd0:levd1,lond0:lond1,latd0:latd1 ) )

  call georeferenced_data_allocate(qset%total)

end subroutine

!> Deallocates the geometric-height arrays, the quantity sets, and the
!! quantities array.
subroutine deallocate_quantity_computation_module()

  implicit none

  if(allocated(zg_int)) deallocate(zg_int)
  if(allocated(zg_mid)) deallocate(zg_mid)
  if(allocated(zg_mid_nm)) deallocate(zg_mid_nm)
  if(allocated(zg_int_nm)) deallocate(zg_int_nm)

  if(allocated(zg_mid_mean)) deallocate(zg_mid_mean)
  if(allocated(zg_mid_mean_nm)) deallocate(zg_mid_mean_nm)
  if(allocated(zg_int_mean)) deallocate(zg_int_mean)
  if(allocated(zg_int_mean_nm))deallocate(zg_int_mean_nm)

  call qset%deallocate
  if(allocated(quantities))deallocate(quantities)

end subroutine

!> Computes geometric height and its ensemble mean, evaluates the requested
!! quantities (full set for analysis steps, output set otherwise), and
!! interpolates them onto all georeferenced datasets.
subroutine compute_results(for_analysis_step)

  ! tie-gcm
  use fields_module, only: itc

  ! intern
  use georeferenced_data_module, only: georeferenced_datasets,n_datasets
  use quantity_info_module, only: STEP_CURRENT, STEP_PREVIOUS, get_info
  use mpi_moments_module, only: compute_mean_mpi

  implicit none

  ! arguments
  logical, intent(in), optional :: for_analysis_step

  ! local
  integer :: i
  logical :: for_analysis_step_

  if(present(for_analysis_step)) then
    for_analysis_step_ = for_analysis_step
  else
    for_analysis_step_ = .false.
  end if


  call compute_geometric_height(itc,zg_mid,zg_int,STEP_CURRENT)
  call compute_geometric_height(itc,zg_mid_nm,zg_int_nm,STEP_PREVIOUS)

  call compute_mean_mpi(zg_int, zg_int_mean)
  call compute_mean_mpi(zg_mid, zg_mid_mean)
  call compute_mean_mpi(zg_int_nm, zg_int_mean_nm)
  call compute_mean_mpi(zg_mid_nm, zg_mid_mean_nm)

  if(allocated(quantities)) deallocate(quantities)

  if(for_analysis_step_)then
    allocate(quantities(qset%total%size()))
    do i=1, qset%total%size()
      quantities(i)%info=>get_info(qset%total%at(i))
    end do
  else
    allocate(quantities(qset%output%size()))
    do i=1, qset%output%size()
      quantities(i)%info=>get_info(qset%output%at(i))
    end do
  end if

  do i=1, size(quantities)
    if(cfg_log%verbose_level>0) write(*,*) 'compute ', quantities(i)%info%name
    call quantities(i)%info%calc()
  end do

  do i=1,n_datasets
    if(cfg_log%verbose_level>0) write(*,*) 'interpolate to ', georeferenced_datasets(i)%ptr%name
    call georeferenced_datasets(i)%ptr%interpolate(quantities, zg_mid, zg_mid_nm, zg_int, zg_int_nm)
  end do

end subroutine

!> Computes geometric height at model levels and half-levels from
!! TN/O2/O1/HE for the given time step (current or previous).
subroutine compute_geometric_height(itx,zg_mid,zg_int,step)

  ! intern
  use aerostatic_diag, only: calc_diag, int_to_mid
  use quantity_info_module, only: STEP_CURRENT, STEP_PREVIOUS

  ! tie-gcm
  use fields_module, only: levd0,levd1,lond0,lond1,latd0,latd1, &
                           f4d, &
                           i_tn, i_o2, i_o1, i_he, i_tn_nm, i_o2_nm, i_o1_nm, i_he_nm

  implicit none

  integer, intent(in) :: itx
  real,  dimension(levd0:levd1,lond0:lond1,latd0:latd1), intent(out) :: &
    zg_mid,zg_int
  integer, intent(in):: step

  ! local
  integer, dimension(4) :: fids

  select case(step)
    case(STEP_CURRENT)
      fids = (/i_tn, i_o2, i_o1, i_he/)
    case(STEP_PREVIOUS)
      fids = (/i_tn_nm, i_o2_nm, i_o1_nm, i_he_nm/)
    case default
      ! set to zero to avoid -Wmaybe-uninitialized warning
      fids = 0
      call shutdown('compute_geometric_height invalid step')
  end select

  call calc_diag(lon0=lond0,lon1=lond1,&
                 lev0=levd0,lev1=levd1,&
                 lat0=latd0,lat1=latd1,&
                 tn=f4d(fids(1))%data(:,:,:,itx), &
                 o2=f4d(fids(2))%data(:,:,:,itx), &
                 o1=f4d(fids(3))%data(:,:,:,itx), &
                 he=f4d(fids(4))%data(:,:,:,itx), &
                 zg=zg_int)

  zg_mid = int_to_mid(lond0,lond1,levd0,levd1,latd0,latd1,zg_int)

end subroutine


end module
