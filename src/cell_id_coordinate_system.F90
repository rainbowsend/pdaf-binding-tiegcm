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
! Converts geographic lon/lat/level coordinates to continuous TIE-GCM grid-cell-index coordinates for interpolation.

module cell_id_coordinate_system
  ! computes the position of observation w.r.t to TIE-GCM grid id starting at (0,0,0)
  !
  !   |--1--|--2--|--3--|--4--|  cell index
  !   0     1     2     3     4  coordinates (cell boundaries)
  !   | 0.5 | 1.5 | 2.5 | 3.5 |  cell center coordinate
  !
  ! for periodic dimension
  !   |--1--|--2--|--3--|--4--|--1--|  cell index
  !   0     1     2     3     4     0  coordinates (cell boundaries)
  !   | 0.5 | 1.5 | 2.5 | 3.5 | 0.5 |  cell center coordinate

  ! tie-gcm
  use fields_module,only: levd0,levd1,lond0,lond1,latd0,latd1

  ! intern
  use interpolation_module, only: spline_interpolator

  implicit none

  real, dimension(:), allocatable, protected :: cell_id_lev ! cell center coordinate in vertical direction
  real, dimension(:), allocatable, protected :: cell_id_lon ! cell center coordinate in zonal direction
  real, dimension(:), allocatable, protected :: cell_id_lat ! cell center coordinate in meridional direction

  type(spline_interpolator), private :: zonal
  type(spline_interpolator), private :: meridional

  interface to_zonal_meridional_vertical
    module procedure to_zonal_meridional_vertical_1, to_zonal_meridional_vertical_2
  end interface

  contains

  !> Initializes the cell-id coordinate arrays and the zonal/meridional spline interpolators used to convert geographic positions to continuous grid-cell coordinates.
  subroutine init_cell_id_coordinate_system

    ! intern
    use tiegcm_optimized_interpolator, only: lon_p_halo, lat_p_halo

    implicit none

    ! local
    integer :: i

    allocate(cell_id_lev(levd0:levd1))
    allocate(cell_id_lon(lond0:lond1))
    allocate(cell_id_lat(latd0:latd1))

    do i=lbound(cell_id_lev,dim=1),ubound(cell_id_lev,dim=1)
      cell_id_lev(i) = i-.5
    end do
    do i=lbound(cell_id_lon,dim=1),ubound(cell_id_lon,dim=1)
      cell_id_lon(i) = i-2.5 ! substract two periodic fields
    end do
    do i=lbound(cell_id_lat,dim=1),ubound(cell_id_lat,dim=1)
      cell_id_lat(i) = i-.5
    end do

    call zonal%init(lon_p_halo,cell_id_lon,degree=1)
    call meridional%init(lat_p_halo,cell_id_lat,degree=1)

  end subroutine

  !> Converts a longitude to its continuous TIE-GCM grid-cell-index coordinate via 1D spline interpolation.
  function longitude_to_cell_id_coordinate(lon) result (cell_id_coord)

    implicit none
    ! arguments
    real, intent(in) :: lon

    ! result
    real :: cell_id_coord

    cell_id_coord = zonal%interpolate(lon)
  end function

  !> Converts a latitude to its continuous TIE-GCM grid-cell-index coordinate via 1D spline interpolation.
  function latitude_to_cell_id_coordinate(lat) result (cell_id_coord)

    implicit none
    ! arguments
    real, intent(in) :: lat

    ! result
    real :: cell_id_coord

    cell_id_coord = meridional%interpolate(lat)
  end function

  !> Converts a geographic (lon, lat, height) point to its (lev, lon, lat) cell-id coordinates, interpolating the vertical index from the given geometric/geopotential height field.
  function to_cell_id_coordinate(zg,zglevel,point) result(cell_id_coordinate)

    use quantity_info_module, only: quantity_info, LEVEL_MID, LEVEL_INT, STEP_CURRENT
    use tiegcm_optimized_interpolator, only: sparse_state_interpolator

    implicit none

    ! arguments
    real, dimension(:,:,:), intent(inout) :: zg ! size as in fields module
    integer, intent(in) :: zglevel ! LEVEL_MID or LEVEL_INT
    real, dimension(3), intent(in) :: point ! lon (-180 deg :180 deg), lat(-90 deg : 90 deg), alt (m)

    ! result
    real, dimension(3) :: cell_id_coordinate

    ! local
    type(sparse_state_interpolator) :: vertical
    real, dimension(levd0:levd1,lond0:lond1,latd0:latd1) :: M
    type(quantity_info) :: info

    cell_id_coordinate = 0.0
    cell_id_coordinate(2) = longitude_to_cell_id_coordinate(point(1))
    cell_id_coordinate(3) = latitude_to_cell_id_coordinate(point(2))

    call cell_id_meshgrid(M,dim=1)

    select case(zglevel)
      case(LEVEL_MID)
        call vertical%init(point, &
                           zg_mid=zg, &
                           degree=1)
      case(LEVEL_INT)
        call vertical%init(point, &
                           zg_int=zg, &
                           degree=1)
    end select

    info%name = "idx"
    info%level=zglevel
    info%step =STEP_CURRENT
    call vertical%interpolate(info, M, cell_id_coordinate(1))
    call vertical%deallocate()

!      write(*,*) ' cell id coord :', cell_id_coordinate, ' geometric coord: ', point

  end function

  !> Fills a 3D array with the cell-id coordinate values broadcast along the given dimension (1=lev, 2=lon, 3=lat).
  subroutine cell_id_meshgrid(M,dim)

    implicit none

    ! arguments
    real, dimension(levd0:levd1,lond0:lond1,latd0:latd1), intent(out) :: M
    integer, intent(in) :: dim

    ! local
    integer :: i

    do i=lbound(M,dim=dim),ubound(M,dim=dim)
      select case (dim)
        case (1) ! lev
          M(i,:,:) = cell_id_lev(i)
        case (2) ! lon
          M(:,i,:) = cell_id_lon(i)
        case (3) ! lat
          M(:,:,i) = cell_id_lat(i)
      end select
    end do
  end subroutine

  !> Reorders a (lev,lon,lat) coordinate triple to (lon,lat,lev).
  subroutine to_zonal_meridional_vertical_1(coords)
    implicit none

    ! arguments
    real, dimension(3), intent(inout) :: coords

    ! local
    real, dimension(3) :: tmp

    tmp=coords

    coords(1) = tmp(2)
    coords(2) = tmp(3)
    coords(3) = tmp(1)

  end subroutine

  !> Reorders each (lev,lon,lat) coordinate triple in a 2D array of coordinates to (lon,lat,lev).
  subroutine to_zonal_meridional_vertical_2(coords)
    implicit none

    ! arguments
    real, dimension(:,:), intent(inout) :: coords

    ! local
    real, dimension(:,:), allocatable :: tmp

    allocate(tmp,source=coords)

    coords(1,:) = tmp(2,:)
    coords(2,:) = tmp(3,:)
    coords(3,:) = tmp(1,:)

    deallocate(tmp)

  end subroutine

end module
