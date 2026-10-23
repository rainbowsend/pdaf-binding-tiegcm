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
! Provides a function for computing vertical electron content from TIE-GCM
! electron densites (NE) considering missing contributiuons from plasmasphere
! and ionosphere above upper TIE-GCM boundary

module vtec_module

  implicit none

  real, allocatable :: vpec(:,:) ! Plasmaspheric vertical electron content

  real, allocatable :: dp_lat(:,:,:) ! magnetic dipole latitude at model grid(rad)
  real, allocatable :: dp_lon(:,:,:) ! magnetic dipole longitude at model grid (rad)
  real, allocatable :: L(:,:,:)    ! apex radius

  integer :: n_vert = 100    ! number of levels for interpolation of plasmasphere
  real, allocatable :: z(:)  ! altitudes of levels where plasmasphere electron number density is computed
  real :: dz                 ! vertical resolution of plasmasphere integration
  real, private :: pZ0       ! lower altitude of selected plasmasphere model

  ! Currently two empirical plasmasphere electron number density models are implemented
  integer, parameter :: PLASMASPHERE_DENSITY_OZHOGIN = 1  ! Ozhogin et al. (2012) https://doi.org/10.1029/2011JA017330
  integer, parameter :: PLASMASPHERE_DENSITY_HARTLEY = 2  ! Hartley et al. (2023) https://doi.org/10.1029/2022JA031012

  integer :: density_model = PLASMASPHERE_DENSITY_HARTLEY

  ! 1 tecu = 1e+16 electrons/m2 = 1e+12 electrons/cm2
  real, parameter :: cm2_per_tecu = 1e+12



contains

subroutine init_vtec

  ! tiegcm
  use cons_module, only: re, dtr
  use params_module, only: glon, glat

  ! intern
  use state_module, only: lonX0, lonX1, latX0, latX1

  use array_print_module, only: printMat

  implicit none

  ! local
  integer :: i

  allocate (dp_lat(lonX0:lonX1, latX0:latX1, n_vert))
  allocate (dp_lon(lonX0:lonX1, latX0:latX1, n_vert))
  allocate (L(lonX0:lonX1, latX0:latX1, n_vert))

  allocate(vpec(lonX0:lonX1,latX0:latX1))

  allocate(z(n_vert))

  select case (density_model)
    case(PLASMASPHERE_DENSITY_OZHOGIN)
      pZ0 = 2000E+5
    case(PLASMASPHERE_DENSITY_HARTLEY)
      pZ0 = 1000E+5
    case default
    call shutdown('unknown plasmaspheric model')
  end select

  dz = 6*re/n_vert
  z(1)=pZ0
  do i=2, n_vert
    z(i) = z(i-1) + dz
  end do

  ! TODO use IGRF derived magnetic coordinates instead of dipole
  call approximate_quasi_dipole_and_apex(glon(lonX0-2:lonX1-2),glat(latX0:latX1),z,dp_lon,dp_lat,L)

!!   For testing
!   call printMat(dp_lon(:,:,1)/dtr,"dp_lon","f7.2")
!   call printMat(dp_lat(:,:,1)/dtr,"dp_lat","f7.2")
!   call printMat(L,"L","f7.2")
!
!   call calc_vpec()
!
!   call printMat(vpec,"vpec")

end subroutine

!> Computes vertical total electron content in tec units
!!
!! Integrates 3D electron number density and extrapolates it up to the lower boundary of the selected plasmasphere model
!! Contribution of plasmasphere is added
function vertical_total_electron_content(ne, zg) result(vtec)

  implicit none

  ! arguments
  real, dimension(:,:,:), intent(in) :: ne ! electron density at interfaces (1/cm3)
  real, dimension(:,:,:), intent(in) :: zg ! geometric height at interfaces (cm)

  ! result
  real, dimension(size(ne,dim=2),size(ne,dim=3)) :: vtec ! vertical electron content (tecu)

  ! integrate electron density and extrapolate vertical height profile up to lowerboundary of plasmasphere model
  vtec = integrate_electron_density(ne, zg, extrapolation_limit=pZ0) / cm2_per_tecu

  ! compute and add contribution from plasmasphere
  call calc_vpec()
  vtec = vtec + vpec

end function vertical_total_electron_content

!> Integrates the electron density over the vertical column.
!!
!! Total electron number within vertical column is approximated via trapezoidal rule.
!! Contribution from ionosphere above upper boundary of model can optionally be approximated
!! assuming exponential decay up to the specified limit for extrapolation
function integrate_electron_density(ne, zg, extrapolation_limit) result(column)

  implicit none

  ! arguments
  real, dimension(:,:,:), intent(in) :: ne ! electron density at interfaces (1/cm3)
  real, dimension(:,:,:), intent(in) :: zg ! geometric height at interfaces (cm)
  real, intent(in) :: extrapolation_limit  ! extrapolate vertical profile beyond upper boundary up to this altitude (cm)
                                           ! if negative do not apply extrapolation

  ! result
  real, dimension(size(ne,dim=2),size(ne,dim=3)) :: column ! electron content (1/cm2)

  ! local
  integer :: nlev
  real, dimension(size(ne,dim=2),size(ne,dim=3)) :: N0       ! electron number density at top of TIE-GCM (1/cm3)
  real, dimension(size(ne,dim=2),size(ne,dim=3)) :: H        ! electron number density scale height (cm)
  real, dimension(size(ne,dim=2),size(ne,dim=3)) :: ne_extra ! Contribution of ionosphere above upper boundary (1/cm2)

  nlev = size(ne,dim=1)
  ! vectorized trapezoidal rule
  ! integral f(x) from a to b = 0.5 (b-a) (f(a)+f(b))
  column = 0.5 * sum(  (zg(2:nlev,:,:)-zg(1:nlev-1,:,:)) * (ne(1:nlev-1,:,:)+ne(2:nlev,:,:)), dim=1)

  if(extrapolation_limit>0) then

    ! assume exponential decay above upper boundary of TIE-GCM
    !
    ! N(Z) = N0 exp(-(Z-Z0)/H)
    !
    ! where Z0 and N0 are the height and electron number density
    ! at the upper boundary
    N0 = ne(nlev,:,:)

    ! Calculate scale height from the two upper most levels
    ! extrapolation with exponential is only meaningful if electron density
    ! decays at the top of the model
    where( ne(nlev-1,:,:) > ne(nlev,:,:) )
      H = (zg(nlev,:,:)-zg(nlev-1,:,:)) / log(ne(nlev-1,:,:)/N0)
    elsewhere
      H = 0
    end where

    ! integral over N(Z) from Z0 to +inf is H N0
    ! ne_extra = H * N0
    ! integral over N(Z) from Z0 to Z1 with Z1>Z0
    ne_extra = 0
    where(zg(nlev,:,:) < extrapolation_limit .and. H > 0)
      ne_extra = N0 * ( H - H*exp( (zg(nlev,:,:)-extrapolation_limit)/H ))
    end where
    column = column + ne_extra

  end if

end function integrate_electron_density

!> Centered-dipole latitude and L-shell of a point
subroutine approximate_quasi_dipole_and_apex(lon,lat,alt,dp_lon,dp_lat,L)

  ! tiegcm
  use cons_module, only: re

  implicit none

  ! arguments
  real, dimension(:), intent(in) :: lon ! geographic lonitdue (deg)
  real, dimension(:), intent(in) :: lat ! geographic latitude (deg)
  real, dimension(:), intent(in) :: alt ! altitude above surface in (cm)
  real, dimension(:,:,:), intent(out) :: dp_lon ! magnetic dipole longitdue (rad)
  real, dimension(:,:,:), intent(out) :: dp_lat ! magnetic dipole latitude (rad)
  real, dimension(:,:,:), intent(out) :: L    ! Apex radius (Earth radii)

  ! local
  integer :: i
  real, dimension(size(lon),size(lat)) :: lons, lats

  lons = spread(lon,2,size(lat))
  lats = spread(lat,1,size(lon))

  call dipole_logitude_and_latitude(lons,lats,dp_lon(:,:,1),dp_lat(:,:,1))

  dp_lon = spread(dp_lon(:,:,1),3,size(alt))
  dp_lat = spread(dp_lat(:,:,1),3,size(alt))

  do i=1, size(alt)
    ! r(dp_lat) = L Re cos^2(dp_lat)
    ! L = r(dp_lat)/(Re cos^2(dp_lat))
    L(:,:,i)  = ((re+alt(i))/re) / cos(dp_lat(:,:,i))**2
  end do

end subroutine

!> First order approximation of Earth's magnetic field assuming dipole structure
! rotate z axis to the magnetic pole
elemental subroutine dipole_logitude_and_latitude(lon,lat,dp_lon,dp_lat)

  ! tiegcm
  use cons_module, only: dtr
  use apex, only: colat, elon

  implicit none

  ! arguments
  real, intent(in) :: lon ! geographic longitude (deg)
  real, intent(in) :: lat ! geographic latitude (deg)
  real, intent(out) :: dp_lon ! magnetic dipole longitude (rad)
  real, intent(out) :: dp_lat ! magnetic dipole latitude (rad)

  real :: lat_pole, lon_pole

  lat_pole = 90. - colat
  lon_pole = elon

  ! Van Zandt et al. 1972 eq. 2
  dp_lat = asin( sin(lat*dtr) * sin(lat_pole*dtr) + cos(lat*dtr) * cos(lat_pole*dtr) * cos(lon*dtr-lon_pole*dtr) )

  ! the equation in the paper maps longitude only from -90 to 90
  !   dp_lon = asin( cos(lat*dtr) * sin((lon-lon_pole)*dtr) / cos(dp_lat) )
  ! using atan2 gives the full circle
  dp_lon = atan2( cos(lat*dtr) * sin((lon-lon_pole)*dtr), &
                  sin(lat_pole*dtr) * cos(lat*dtr)*cos((lon-lon_pole)*dtr) - cos(lat_pole*dtr)*sin(lat*dtr) )


end subroutine

subroutine calc_vpec()

  ! tiegcm
  use cons_module, only: dtr
  use magfield_module, only: sunlons

  ! intern
  use state_module, only: lonX0, lonX1, latX0, latX1
!   use nc_debug, only: instant_write_2d

  implicit none

  ! local
  real, dimension(lonX0:lonX1, latX0:latX1, n_vert) :: ne ! electron number density (1/cm3)
  real, dimension(lonX0:lonX1, latX0:latX1) :: mlt ! magnetic local time (hours)

  ne = 0

   select case (density_model)
    case(PLASMASPHERE_DENSITY_OZHOGIN)
      where(1.5 < L .and. L < 4.)
        ne = ozhogin_density(L,dp_lat)
      endwhere
    case(PLASMASPHERE_DENSITY_HARTLEY)
      mlt = modulo(12. + (dp_lon(:,:,1) - sunlons(1))/dtr/15., 24.)

      where(1. < L .and. L < 3.)
        ne = hartley_density(L,dp_lat,spread(mlt,3,n_vert))
      endwhere
    case default
      call shutdown('unknown plasmaspheric model')
  end select

  vpec = dz*sum(ne,dim=3)/cm2_per_tecu

!   call instant_write_2d("vpec", vpec)

end subroutine


!> Plasmaspheric electron density (1/cm3), Hartley et al. (2023) eq. 6-8.
!! https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2022JA031012
elemental function hartley_density(L,mlat,mlt) result (ne)

  ! tiegcm
  use cons_module, only: pi

  implicit none

  ! arguments
  real, intent(in) :: L    ! apex radius, shell parameter (Earth radii)
  real, intent(in) :: mlat ! magnetic dipole latitude (rad)
  real, intent(in) :: mlt ! magnetic local time (hours)

  ! result
  real :: ne ! (1/cm3)

  ! local
  real :: ne_eq
  real :: A
  real :: alat ! apex latitude (rad)

  alat = sign(acos(1/sqrt(L)), mlat)

  if (L>1.45) then
    ne_eq = 10**(4.4693-0.4903*L)
  else
    ne_eq = 10**(4.03-2.05*(L-1.24)+3.76*(L-1.24)**2)
  end if

  A = 0.15 + 0.6 * exp(-0.5 * (5*(L-1))**2 )

  ne = ne_eq * ( cos(pi/2 * 1.03*mlat/alat) )**(-0.44) &
             * ( 1 + A *sin(2*pi * (mlt-8)/24) )

end function

!> Plasmaspheric electron density (1/cm3), Ozhogin et al. (2012) eq. 9.
elemental function ozhogin_density(L,mlat) result (ne)

  ! tiegcm
  use cons_module, only: pi

  implicit none

  ! arguments
  real, intent(in) :: L    ! apex radius, shell parameter (Earth radii)
  real, intent(in) :: mlat ! magnetic latitude (rad)

  ! result
  real :: ne ! (1/cm3)

  ! local
  real :: alat ! apex latitude (rad)

  alat = sign(acos(1/sqrt(L)), mlat)

  ! Equation 9 in P. Ozhogin (2012) https://doi.org/10.1029/2011JA017330
  ne = 10**(4.4693-0.4903*L) * cos(pi/2 * mlat/alat)**(-0.75)

end function

end module
