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
! Generic regular-grid subdomain partitioning/index-mapping utility used for localized analysis domains.

module structured_gird_subdomain_module
  implicit none

  integer, protected :: n_gird_dim1 ! size first dimension
  integer, protected :: n_gird_dim2 ! size second dimension
  integer, protected :: n_gird_dim3 ! size third dimension

  integer, protected :: nsub1 ! number subdomains in first dim
  integer, protected :: nsub2 ! number subdomains in second dim
  integer, protected :: nsub3 ! number subdomains in third dim

  integer, protected  :: nsub ! number of subdomains

  integer, protected, allocatable, dimension(:) :: count_1
  integer, protected, allocatable, dimension(:) :: count_2
  integer, protected, allocatable, dimension(:) :: count_3

  integer, protected, allocatable, dimension(:) :: offset_1
  integer, protected, allocatable, dimension(:) :: offset_2
  integer, protected, allocatable, dimension(:) :: offset_3

  type subdomain
    integer :: id
    integer, allocatable, dimension(:) :: mapping  ! indices of state vector used in subdomain
    integer :: domain_shape(3)
    integer :: domain_size
    integer :: offset(3)
    integer :: center_id(3)  ! starting at (1,1,1)
  end type subdomain

  type(subdomain), protected, allocatable, dimension(:) :: subdomains


! if coresponding count (size of subdomain) is uneven the mapping index is exactly the center
! else the id is the cell left of the center, so the the value of the center can be determined easily with linear interpolation using id and id+1
!
!     | | |X| | |         | |X| | |         | |X| |           |X| |          |X|


  contains

  !> Deallocates a single subdomain's index-mapping array.
  subroutine deallocate_subdomain(this)

    implicit none
    type(subdomain), intent(inout) :: this

    if (allocated (this%mapping)) deallocate (this%mapping)
  end subroutine

  !> Prints a subdomain's shape, size, offset, and center-cell index.
  subroutine print_subdomain(this)

    implicit none
    type(subdomain), intent(in) :: this

    write(*,"(i5,'|',4i5,'|',3i5,'|',3i5)") &
      this%id,this%domain_shape, this%domain_size, this%offset, this%center_id

  end subroutine

  !> Deallocates all module-level subdomain partitioning arrays and the subdomains
  !! array itself.
  subroutine deallocate_structured_gird_subdomain_module()

    implicit none

    integer :: i

    if (allocated (count_1)) deallocate (count_1)
    if (allocated (count_2)) deallocate (count_2)
    if (allocated (count_3)) deallocate (count_3)
    if (allocated (offset_1)) deallocate (offset_1)
    if (allocated (offset_2)) deallocate (offset_2)
    if (allocated (offset_3)) deallocate (offset_3)

    if (allocated (subdomains)) then
      do i=1,nsub
        call deallocate_subdomain(subdomains(i))
      end do
      deallocate (subdomains)
    end if

  end subroutine

  !> Sets the full (undecomposed) grid size in each of the three dimensions.
  subroutine init_structured_gird_subdomain_module(n_1,n_2,n_3)
    integer, intent(in) :: n_1 ! size
    integer, intent(in) :: n_2 ! size
    integer, intent(in) :: n_3 ! size

    n_gird_dim1=n_1
    n_gird_dim2=n_2
    n_gird_dim3=n_3
  end subroutine


  !> Partitions the grid into disjoint subdomains of approximately size (d1,d2,d3)
  !! in each dimension, computing per-subdomain counts and offsets.
  subroutine disjoint_domain_regulary(d1,d2,d3)

    use task_distribution_module, only: get_offset, distribute

    integer, intent(in) :: d1 ! size of subdomains in first dim
    integer, intent(in) :: d2 ! size of subdomains in second dim
    integer, intent(in) :: d3 ! size of subdomains in third dim

    nsub1 = CEILING(real(n_gird_dim1)/real(d1))
    nsub2 = CEILING(real(n_gird_dim2)/real(d2))
    nsub3 = CEILING(real(n_gird_dim3)/real(d3))

    allocate(offset_1(nsub1))
    allocate(offset_2(nsub2))
    allocate(offset_3(nsub3))
    allocate(count_1(nsub1))
    allocate(count_2(nsub2))
    allocate(count_3(nsub3))

    call distribute( nsub1, n_gird_dim1, count_1 )
    call distribute( nsub2, n_gird_dim2, count_2 )
    call distribute( nsub3, n_gird_dim3, count_3 )

    offset_1 = get_offset( count_1, 1 )
    offset_2 = get_offset( count_2, 1 )
    offset_3 = get_offset( count_3, 1 )


  end subroutine disjoint_domain_regulary

  !> Builds the index mapping from each subdomain to its flattened grid-cell
  !! indices, given the previously computed subdomain counts/offsets.
  subroutine compute_mapping()
    ! computes the mapping given counts and offsets

    ! intern
    use array_print_module, only: printMat

    integer, allocatable, dimension(:,:,:) :: idx
    integer, allocatable, dimension(:) :: idx_flat

    integer:: i,j,k,c
    integer:: n

    integer, allocatable, dimension(:,:,:) :: sub_idx

    nsub=nsub1*nsub2*nsub3

    allocate( subdomains(nsub) )

    n=n_gird_dim1*n_gird_dim2*n_gird_dim3

    allocate(sub_idx(maxval(count_1),maxval(count_2),maxval(count_3)))

    allocate(idx_flat(n))
    do i=1,n
      idx_flat(i)=i
    end do

    allocate(idx(n_gird_dim1, n_gird_dim2, n_gird_dim3))
    idx(:,:,:) = reshape( idx_flat, (/n_gird_dim1, n_gird_dim2, n_gird_dim3 /) )
!     call printMat(idx)


     write(*,"(a)") "subdomain indexing"
     write(*,"(a5,'|',4a5,'|',3a5,'|',3a5)") &
      "id", "cnt1", "cnt2", "cnt3", "size",  "off1", "off2", "off3", "cid1", "cid2", "cid3" 

    c=1
    do i = 1,nsub1
      do j = 1,nsub2
        do k = 1,nsub3
          subdomains(c)%id = c
          subdomains(c)%domain_size = count_1(i)*count_2(j)*count_3(k)
          allocate(subdomains(c)%mapping( subdomains(c)%domain_size ))

          sub_idx(1:count_1(i),1:count_2(j),1:count_3(k)) = &
            idx(offset_1(i):offset_1(i)+count_1(i)-1,&
                offset_2(j):offset_2(j)+count_2(j)-1,&
                offset_3(k):offset_3(k)+count_3(k)-1)

!           call printMat(sub_idx(1:count_1(i),1:count_2(j),1:count_3(k)))

          subdomains(c)%mapping(:) = &
            reshape( sub_idx(1:count_1(i),1:count_2(j),1:count_3(k)),&
            (/subdomains(c)%domain_size/) )

          ! index of center cell
          subdomains(c)%center_id(1) = offset_1(i)+ceiling(count_1(i)/2.)-1
          subdomains(c)%center_id(2) = offset_2(j)+ceiling(count_2(j)/2.)-1
          subdomains(c)%center_id(3) = offset_3(k)+ceiling(count_3(k)/2.)-1

          subdomains(c)%domain_shape(1) = count_1(i)
          subdomains(c)%domain_shape(2) = count_2(j)
          subdomains(c)%domain_shape(3) = count_3(k)

          subdomains(c)%offset(1) = offset_1(i)
          subdomains(c)%offset(2) = offset_2(j)
          subdomains(c)%offset(3) = offset_3(k)

          call print_subdomain(subdomains(c))

          c=c+1
        end do
      end do
    end do

    deallocate(idx_flat)
    deallocate(sub_idx)
    deallocate(idx)

  end subroutine compute_mapping

  !> Scatters a subdomain-local array of values into their corresponding positions
  !! in the full local state vector.
  subroutine insert_local_subdomain(domain_id,val,state_p)

    implicit none

    ! arguments
    integer, intent(in) :: domain_id
    real,dimension(:),intent(in) :: val
    real,dimension(:),intent(inout) :: state_p

    ! local
    integer :: i

    do i=1,subdomains(domain_id)%domain_size
      state_p(subdomains(domain_id)%mapping(i))=val(i)
    end do
  end subroutine

  !> Gathers a subdomain's values from the full local state vector into a
  !! subdomain-local array.
  subroutine extract_local_subdomain(domain_id,val,state_p)

    implicit none

    ! arguments
    integer, intent(in) :: domain_id
    real,dimension(:),intent(inout) :: val
    real,dimension(:),intent(in) :: state_p

    ! local
    integer :: i

    do i=1,subdomains(domain_id)%domain_size
      val(i)=state_p(subdomains(domain_id)%mapping(i))
    end do
  end subroutine

  !> Test/debug routine: computes and prints the center coordinates of every
  !! subdomain.
  subroutine test_get_center_coords

    ! intern
    use array_print_module, only: printMat

    implicit none

    ! local
    integer :: i,j,k,c
    real, dimension(nsub1,nsub2,nsub3,3) :: center

    ! same order as in compute_mapping
    c=1
    do i = 1,nsub1
      do j = 1,nsub2
        do k = 1,nsub3
          center(i,j,k,:) = get_center_coords_cell_id(c,0,0)
          c=c+1
        end do
      end do
    end do

    call printMat(center(:,:,:,1),name='lev',format='f4.0')
    call printMat(center(:,:,:,2),name='lon',format='f4.0')
    call printMat(center(:,:,:,3),name='lat',format='f4.0')

  end subroutine


  !> Interpolates the geographic center coordinates (lon, lat, altitude) of a
  !! subdomain's center cell from a TIE-GCM coordinate field Z. TIE-GCM specific.
  function get_center_coords(domain_id,lon,lat,Z) result(center)
  ! specialised function for TIEGCM
  ! center(1) longitude in radian
  ! center(2) latitude  in radian
  ! center(3) altitude  in meters

    ! intern
    use array_print_module, only: printMat
    use math_addon_module, only: pi

    implicit none

    integer, intent(in) :: domain_id
    real,dimension(:),intent(in) :: lon, lat
    real,dimension(:,:,:),intent(in) :: Z

    real,dimension(3) :: center ! center coordinates of subdomain (lon,lat,alt)

    ! local
    real,dimension(2,2) :: interp_lon
    real,dimension(2) :: interp_lat

    integer :: c
    integer :: center_id(3)
    integer :: lev_pair(2) ! central grid point and its upper neighbour
    integer :: lat_pair(2) ! central grid point and its upper neighbour

    ! alias
    c=domain_id
    center_id=subdomains(c)%center_id

    ! The upper neighbour is required for an even extent only, but the slices
    ! below are read in either case. It therefore has to stay inside the array:
    ! a subdomain of extent one at the upper boundary of the grid does not have
    ! one. Clamping does not change the result, since for an even extent the
    ! neighbour always exists and for an odd extent its value is discarded.
    lev_pair = (/ center_id(1), min(center_id(1)+1, size(Z,dim=1)) /)
    lat_pair = (/ center_id(3), min(center_id(3)+1, size(Z,dim=3)) /)

    ! order for TIEGCM lev, lon, lat
!
!     call printMat(Z(lev_pair, center_id(2):center_id(2)+1, lat_pair), "Z")


    ! get longitude of the center of the subdomain.
    ! get latitude - level slice of Z at central longitude
    if( mod(subdomains(c)%domain_shape(2),2)==0 )then
      ! even zonal extend of subdomain
      ! compute the longitude between both central grid points
      center(1)= (lon(center_id(2)) + lon(center_id(2)+1))/2
      interp_lon = (Z(lev_pair, center_id(2)  , lat_pair)&
                   +Z(lev_pair, center_id(2)+1, lat_pair))/2
    else
      ! odd zonal extend of subdomain:
      ! just take the coordinate of the central grid point in the subdomain
      center(1)=lon(center_id(2))
      interp_lon= Z(lev_pair, center_id(2), lat_pair)
    end if

    ! interp_lon: lev x lat
!     call printMat(interp_lon,"lon")

    ! get latitude of the center of the subdomain.
    ! get level slice of Z at central longitude and latitude
    if( mod(subdomains(c)%domain_shape(3),2)==0 )then
      ! even meridional extend of subdomain
      center(2)= (lat(center_id(3)) + lat(center_id(3)+1))/2
      interp_lat=(interp_lon(:,1)+interp_lon(:,2))/2
    else
       ! odd meridional extend of subdomain
      center(2)=lat(center_id(3))
      interp_lat=interp_lon(:,1)
    end if

    ! compute geometric height of subdomain by interpolating over altitude
    if( mod(subdomains(c)%domain_shape(1),2)==0 )then
      ! even vertical extend of subdomain
      center(3)= (interp_lat(1)+interp_lat(2))/2
    else
       ! odd vertical extend of subdomain
      center(3)=interp_lat(1)
    end if

!     write(*,*) center(3),  center(1), center(2)
!     call print_subdomain(subdomains(c))
!     write(*,*) "------------------------------"

    center(1)=center(1)/180.*pi  ! degree to radian
    center(2)=center(2)/180.*pi  ! degree to radian
    center(3)= center(3)/100.    ! cm to meters

  end function

 !> Returns a subdomain's center-cell index coordinates (with optional lon/lat
 !! offsets applied), in the TIE-GCM PDAF cell-index coordinate system. TIE-GCM
 !! specific.
 function get_center_coords_cell_id(domain_id, offset_lon, offset_lat) result(center)

    ! intern
    use array_print_module, only: printMat

    implicit none

    ! arguments
    integer, intent(in) :: domain_id
    integer, intent(in) :: offset_lon
    integer, intent(in) :: offset_lat

    ! result
    real,dimension(3) :: center ! tuple (lev,lon,lat)

    ! local


    center = subdomains(domain_id)%center_id

    ! add offsets due to lon,lat subdivided model grid (parallelization)
    center(2) = center(2)+offset_lon
    center(3) = center(3)+offset_lat

    ! For an even number of cells in a subdomain dimension the center index corresponds to the index left from the real center
    !
    ! The subdomain center coordinate system starts at 1. The coordinate indicates the center of the cell
    ! For the TIE-GCM PDAF system the cell coordinate system starts at 0 and refers to the lower boundary of the cell:
    !
    !   |--1--|--2--|--3--|--4--|  cell index
    !   0     1     2     3     4  coordinates (cell boundaries)
    !   | 0.5 | 1.5 | 2.5 | 3.5 |  cell center coordinate
    !
    ! For an even number of cells in a subdomain dimension the center variable holds the correct information
    ! For an odd number of cells in a subdomain dimension we need to substract -.5
    !
    if( mod(subdomains(domain_id)%domain_shape(1),2)==1 ) center(1)=center(1)-.5
    if( mod(subdomains(domain_id)%domain_shape(2),2)==1 ) center(2)=center(2)-.5
    if( mod(subdomains(domain_id)%domain_shape(3),2)==1 ) center(3)=center(3)-.5

  end function

end module structured_gird_subdomain_module
