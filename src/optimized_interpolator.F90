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
! Vectorized interpolators from TIE-GCM's structured grid to regular/point grids, used to map model state to observation space.

MODULE tiegcm_optimized_interpolator
  ! This module is optimized to interpolate the TIE-GCM structured grid
  ! to a regular grid
  ! first interpolate in veritcal (irregular) dimension
  ! Interpolation in remaining dimensions is vectorized

  ! tiegcm
  use fields_module, only: nf4d

  ! intern
  use bspline_module, only: bspline, sparse_band_matrix
  use configuration, only: cfg_log

  implicit none

  type :: reg_grid_spline_interpolator
    real, dimension(:), allocatable, private :: alt_dst
    real, dimension(:), allocatable, private :: lon_dst
    real, dimension(:), allocatable, private :: lat_dst

    integer :: degree=3

    type(bspline), dimension(:,:), allocatable :: alt_spline
    type(sparse_band_matrix), dimension(:,:), allocatable :: alt_jacobian
    type(bspline) :: lon_spline, lat_spline
    type(sparse_band_matrix) :: lon_jacobian, lat_jacobian

    integer :: n_lon_src
    integer :: n_lat_src

    contains
    procedure, pass(this) :: init => reg_grid_spline_interpolator_init
    procedure, pass(this) :: destroy => reg_grid_spline_interpolator_destroy
    procedure, pass(this) :: set_geometric_height => reg_grid_spline_interpolator_set_geometric_height
    procedure, pass(this) :: interpolate => reg_grid_spline_interpolator_interpolate
  end type

  ! Same as reg_grid_spline_interpolator, but for quantities with no vertical
  ! extent (e.g. VTEC): no alt_dst/alt_spline/alt_jacobian, and no
  ! set_geometric_height, since there is nothing vertical to configure.
  type :: reg_grid_spline_horizontal_interpolator
    real, dimension(:), allocatable, private :: lon_dst
    real, dimension(:), allocatable, private :: lat_dst

    integer :: degree=3

    type(bspline) :: lon_spline, lat_spline
    type(sparse_band_matrix) :: lon_jacobian, lat_jacobian

    integer :: n_lon_src
    integer :: n_lat_src

    contains
    procedure, pass(this) :: init => reg_grid_spline_horizontal_interpolator_init
    procedure, pass(this) :: destroy => reg_grid_spline_horizontal_interpolator_destroy
    procedure, pass(this) :: interpolate => reg_grid_spline_horizontal_interpolator_interpolate
  end type

  type :: sparse_interpolator
    integer :: degree=3
    integer :: lon_first, lon_last
    integer :: lat_first, lat_last
    logical :: on_this_rank
    integer :: subdomain_rank
    real :: lon_dst
    real :: lat_dst
    real :: alt_dst
    type(bspline), dimension(:,:), allocatable :: alt_spline
    type(sparse_band_matrix), dimension(:,:), allocatable :: alt_jacobian
    type(bspline) :: lon_spline, lat_spline
    type(sparse_band_matrix) :: lon_jacobian, lat_jacobian
    contains
    procedure, pass(this) :: init => sparse_interpolator_init
    procedure, pass(this) :: destroy => sparse_interpolator_destroy
    procedure, pass(this) :: set_geometric_height => sparse_interpolator_set_geometric_height
    procedure, pass(this) :: interpolate => sparse_interpolator_interpolate
  end type

  real, allocatable, dimension(:), protected :: lon_p_halo, lat_p_halo ! longitudes and latitudes on current subdomain including halos

  real, dimension(:), allocatable :: lb_lon, ub_lon
  real, dimension(:), allocatable :: lb_lat, ub_lat

  logical, protected :: at_north_end, at_south_end
  logical, protected :: at_west_end, at_east_end

! TODO add polymorphism to get rid of code dublication
!   type, abstract :: state_interpolator_interface
!   end type

  type :: state_interpolator
    type(reg_grid_spline_interpolator) :: inter_int
    type(reg_grid_spline_interpolator) :: inter_int_nm
    type(reg_grid_spline_interpolator) :: inter_mid
    type(reg_grid_spline_interpolator) :: inter_mid_nm
    ! Quantities without vertical extent (LEVEL_NONE) do not depend on the
    ! geometric height. Thus a single interpolator is sufficient. It is
    ! neither specific to mid/int points nor to the current/previous step.
    type(reg_grid_spline_horizontal_interpolator) :: inter_horizontal
    logical, dimension(nf4d) :: prognostic_field_mask ! true if prognostic field was interpolated (set halo to zero afterwards)
    contains
    procedure, pass(this) :: init => state_interpolator_init
    procedure, pass(this) :: interpolate => state_interpolator_interpolate
    procedure, pass(this) :: deallocate => state_interpolator_destroy
  end type state_interpolator

  type :: sparse_state_interpolator
    type(sparse_interpolator) :: inter_int
    type(sparse_interpolator) :: inter_int_nm
    type(sparse_interpolator) :: inter_mid
    type(sparse_interpolator) :: inter_mid_nm
    logical, dimension(nf4d) :: prognostic_field_mask ! true if prognostic field was interpolated (set halo to zero afterwards)
    logical :: on_this_rank
    logical :: exchange_halo
    contains
    procedure, pass(this) :: init => sparse_state_interpolator_init
    procedure, pass(this) :: interpolate => sparse_state_interpolator_interpolate
    procedure, pass(this) :: deallocate => sparse_state_interpolator_destroy
  end type sparse_state_interpolator

  interface is_within_sub_domain
    procedure is_within_sub_domain_scalar, &
              is_within_sub_domain_array
  end interface is_within_sub_domain


  CONTAINS

  !> Initializes per-rank longitude/latitude halo coordinate arrays and
  !! gathers each subdomain's lon/lat bounding box across all MPI ranks.
  subroutine init_tiegcm_optimized_interpolator

    ! extern
    use mpi_f08

    ! tie-gcm
    use fields_module, only: lond0,lond1,latd0,latd1
    use params_module, only: dlon, dlat
    use mpi_module, only:  mytidj, mytidi, mytid, ntaskj, ntaski, TIEGCM_WORLD
    use mod_parallel_pdaf, only: local_ntask

    ! intern
    use state_module,only: latX0, latX1, lonX0, lonX1, lon_p, lat_p

    implicit none

    integer :: i

    write(*,*) 'initalizing interpolation module'

    ! TODO deallocation at end

    allocate(lon_p_halo(lond0:lond1))
    allocate(lat_p_halo(latd0:latd1))

    lon_p_halo(lonX0:lonX1) = lon_p
    do i = lonX0-1,lond0,-1
      lon_p_halo(i) = lon_p_halo(i+1)-dlon
    end do
    do i = lonX1+1,lond1
      lon_p_halo(i) = lon_p_halo(i-1)+dlon
    end do

    lat_p_halo(latX0:latX1) = lat_p
    do i = latX0-1,latd0,-1
      lat_p_halo(i) = lat_p_halo(i+1)-dlat
    end do
    do i = latX1+1,latd1
      lat_p_halo(i) = lat_p_halo(i-1)+dlat
    end do

    write(*,'(a,*(f7.1))') 'lon_p_halo', lon_p_halo
    write(*,'(a,*(f6.1))') 'lat_p_halo', lat_p_halo

    allocate(lb_lon(0:local_ntask-1))
    allocate(ub_lon(0:local_ntask-1))
    allocate(lb_lat(0:local_ntask-1))
    allocate(ub_lat(0:local_ntask-1))

    at_north_end = (mytidj == ntaskj-1)
    at_south_end = (mytidj == 0)

    at_west_end = (mytidi == 0)
    at_east_end = (mytidi == ntaski-1)

    write(*,*) 'at north end ', at_north_end
    write(*,*) 'at south end ', at_south_end
    write(*,*) 'at west end ', at_west_end
    write(*,*) 'at east end ', at_east_end

    ! in case of one PE per model it can contain both north and south end
    lb_lat(mytid) = lat_p_halo(latX0)
    ub_lat(mytid) = lat_p_halo(latX1+1)
    if(at_north_end) ub_lat(mytid) = 90.
    if(at_south_end) lb_lat(mytid) = -90.

    lb_lon(mytid) = lon_p_halo(lonX0)
    ub_lon(mytid) = lon_p_halo(lonX1+1)

    write(*,'(a,f7.1,a,f7.1)') 'range for lons', lb_lon(mytid), '< lon <=',ub_lon(mytid)
    write(*,'(a,f7.1,a,f7.1)') 'range for lats', lb_lat(mytid), '< lat <=',ub_lat(mytid)

    call MPI_Allgather(lb_lon(mytid),1,MPI_REAL8,&
                       lb_lon,1,MPI_REAL8,&
                       TIEGCM_WORLD)
    call MPI_Allgather(ub_lon(mytid),1,MPI_REAL8,&
                       ub_lon,1,MPI_REAL8,&
                       TIEGCM_WORLD)
    call MPI_Allgather(lb_lat(mytid),1,MPI_REAL8,&
                       lb_lat,1,MPI_REAL8,&
                       TIEGCM_WORLD)
    call MPI_Allgather(ub_lat(mytid),1,MPI_REAL8,&
                       ub_lat,1,MPI_REAL8,&
                       TIEGCM_WORLD)

!     do i=0,local_ntask-1
!       write(*,*) 'rank ', i, ' range for lons', lb_lon(i), '< lon <=',ub_lon(i)
!       write(*,*) 'rank ', i, ' range for lats', lb_lat(i), '< lat <=',ub_lat(i)
!     end do


  end subroutine

  !> Initializes the mid/int, current/previous spline interpolators of a
  !! state_interpolator for the given ESMF destination grid, setting
  !! geometric height for each requested component.
  subroutine state_interpolator_init(this, dstField,  zg_mid, zg_mid_nm, zg_int, zg_int_nm, init_flags, degree, grid )

    ! extern
    use ESMF

    ! tie-gcm
    use fields_module,only: levd0,levd1,lond0,lond1,latd0,latd1

    ! intern
    use state_module,only: levX1

    implicit none

    ! arguments
    class(state_interpolator) :: this
    type(ESMF_Field), intent(in), optional :: dstField
    real, dimension(levd0:levd1,lond0:lond1,latd0:latd1), intent(inout), optional :: &
          zg_mid, zg_mid_nm, zg_int, zg_int_nm
    logical, intent(in), optional :: init_flags(4)
    integer, intent(in), optional :: degree
    type(ESMF_Grid), intent(in), optional :: grid

    ! local
    integer :: rc
    type(ESMF_Grid) :: grid_
    real(ESMF_KIND_R8), pointer :: lon_dst(:), lat_dst(:), alt_dst(:)

    logical :: init_flags_(4)

    if( present(init_flags) .eqv. .false.) then
      init_flags_ = .true.
    else
      init_flags_ = init_flags
    end if

    if(present(grid))then
      grid_ = grid
    else
      call ESMF_FieldGet(dstField, grid=grid_, rc=rc)
    end if

    call ESMF_GridGetCoord(grid_,coordDim=1,farrayPtr=lon_dst,rc=rc)
    call ESMF_GridGetCoord(grid_,coordDim=2,farrayPtr=lat_dst,rc=rc)
    call ESMF_GridGetCoord(grid_,coordDim=3,farrayPtr=alt_dst,rc=rc)

    ! ATTENTION upper TIE-GCM layer contains invalid values. Thus we use 1:levX1

    if( (present(zg_mid).eqv..true.) .and. (init_flags_(1).eqv..true.) )then
      if(cfg_log%verbose_level>0) write(*,*) 'init current midlevel points'
      call this%inter_mid%init(lon_p_halo,lat_p_halo, lon_dst, lat_dst, alt_dst, degree)
      call fill_halo_and_rim(zg_mid)
      call this%inter_mid%set_geometric_height(zg_mid(1:levX1,:,:))
    end if
    if( (present(zg_mid_nm).eqv..true.) .and. (init_flags_(2).eqv..true.) )then
      if(cfg_log%verbose_level>0) write(*,*) 'init previous midlevel points'
      call this%inter_mid_nm%init(lon_p_halo,lat_p_halo, lon_dst, lat_dst, alt_dst, degree)
      call fill_halo_and_rim(zg_mid_nm)
      call this%inter_mid_nm%set_geometric_height(zg_mid_nm(1:levX1,:,:))
    end if
    if( (present(zg_int).eqv..true.) .and. (init_flags_(3).eqv..true.) )then
      if(cfg_log%verbose_level>0) write(*,*) 'init current interface points'
      call this%inter_int%init(lon_p_halo,lat_p_halo, lon_dst, lat_dst, alt_dst, degree)
      call fill_halo_and_rim(zg_int)
      call this%inter_int%set_geometric_height(zg_int(1:levX1,:,:))
    end if
    if( (present(zg_int_nm).eqv..true.) .and. (init_flags_(4).eqv..true.) )then
      if(cfg_log%verbose_level>0) write(*,*) 'init previous interface points'
      call this%inter_int_nm%init(lon_p_halo,lat_p_halo, lon_dst, lat_dst, alt_dst, degree)
      call fill_halo_and_rim(zg_int_nm)
      call this%inter_int_nm%set_geometric_height(zg_int_nm(1:levX1,:,:))
    end if

    ! Quantities without vertical extent do not require a geometric height.
    ! Thus this interpolator is always initalized. It is cheap compared to the
    ! ones above, which hold one vertical spline per source column.
    call this%inter_horizontal%init(lon_p_halo,lat_p_halo, lon_dst, lat_dst, degree)

    this%prognostic_field_mask = .false.

  end subroutine

  !> Fills the halo of src_field, marks it if prognostic, and interpolates
  !! it onto dst_field using the matching mid/int, current/previous
  !! component.
  subroutine state_interpolator_interpolate(this, info, src_field, dst_field, logz)

    ! extern
    use ESMF

    ! tie-gcm
    use fields_module,only: levd0,lond0,latd0, f4d

    ! intern
    use array_mapping_module, only: map_to_3d
    use array_print_module, only: printMat
    use character_routines_module, only: string_ends_with
    use quantity_info_module
    use state_module, only: levX1

    implicit none

    ! arguments
    class(state_interpolator) :: this
    type(quantity_info) :: info
    ! assumed shape in the vertical, since quantities without vertical extent
    ! (LEVEL_NONE) are stored with a single level only
    real, dimension(levd0:,lond0:,latd0:), intent(inout) :: src_field
    type( ESMF_Field ), intent(in) :: dst_field
    logical, intent(in), optional :: logz

    ! local
    real(ESMF_KIND_R8), dimension(:,:,:), contiguous, pointer :: dst_ptr ! (lon,lat,alt)

    integer :: rc

    logical :: logz_

    integer :: idx

    if(present(logz)) then
      logz_=logz
    else
      logz_ = .false.
      select case (trim(info%name))
        case ("DEN", "DEN_NM")
          logz_ = .true.
      end select
    end if

    ! the vertical extent of the source field and the level of the quantity
    ! have to be consistent
    if((info%level==LEVEL_NONE) .neqv. (size(src_field,dim=1)==1)) then
      call shutdown('state_interpolator_interpolate: vertical extent of the source field'// &
                    ' does not match the level of quantity '//trim(info%name))
    end if

    ! ATTENTION after interpolation, halo of prognostic fields should be set back to zero,
    ! else the state is slighly changed in the next step. This is done in destoy subroutine
    if(info%level == LEVEL_NONE) then
      ! Exchanging the levels of a quantity without vertical extent would only
      ! be a waste of communication.
      call fill_halo_and_rim_2d(src_field(levd0,:,:))
    else
      call fill_halo_and_rim(src_field)
    end if

    idx = findloc(f4d%short_name,info%name,dim=1)
    if(idx>0) then
      this%prognostic_field_mask(idx)=.true.
    end if

    call ESMF_FieldGet(field=dst_field, localDe=0, farrayPtr=dst_ptr, rc=rc)

    select case(info%level)
      case(LEVEL_MID)
        select case(info%step)
          case(STEP_CURRENT)
            call this%inter_mid%interpolate(src_field(levd0:levX1,:,:),dst_ptr,logz_)
          case(STEP_PREVIOUS)
            call this%inter_mid_nm%interpolate(src_field(levd0:levX1,:,:),dst_ptr,logz_)
        end select
      case(LEVEL_INT)
        select case(info%step)
          case(STEP_CURRENT)
            call this%inter_int%interpolate(src_field(levd0:levX1,:,:),dst_ptr,logz_)
          case(STEP_PREVIOUS)
            call this%inter_int_nm%interpolate(src_field(levd0:levX1,:,:),dst_ptr,logz_)
        end select
      case(LEVEL_NONE)
        ! no vertical extent, thus no distinction between mid/int points and
        ! between current/previous step is required
        call this%inter_horizontal%interpolate(src_field(levd0:levd0,:,:),dst_ptr)
      case default
        call shutdown('state_interpolator_interpolate: unhandled level of quantity '//trim(info%name))
    end select

!     call printMat(dst_ptr,'interp',order=(/3,2,1/))

  end subroutine

  !> Destroys the interpolator components and zeroes the halo of any
  !! prognostic field that was interpolated.
  subroutine state_interpolator_destroy(this)

    ! tie-gcm
    use fields_module, only: f4d,itc

    implicit none

    ! arguments
    class(state_interpolator) :: this

    ! local
    integer :: fidx

    call this%inter_int%destroy
    call this%inter_int_nm%destroy
    call this%inter_mid%destroy
    call this%inter_mid_nm%destroy
    call this%inter_horizontal%destroy

    do fidx = 1, size(this%prognostic_field_mask)
      if(this%prognostic_field_mask(fidx).eqv..true.)then
        if(cfg_log%verbose_level>0) write(*,*) f4d(fidx)%short_name, " set halo zero"
        call zero_halo_and_rim(f4d(fidx)%data(:,:,:,itc))
      end if
    end do

  end subroutine

  !> Initializes a sparse_state_interpolator for a single (lon,lat,alt)
  !! point: determines the owning subdomain rank, whether halo exchange is
  !! needed, and sets up the requested mid/int, current/previous sparse
  !! interpolators.
  subroutine sparse_state_interpolator_init(this, point,  zg_mid, zg_mid_nm, zg_int, zg_int_nm, init_flags, degree )

    ! extern
    use mpi_f08

    ! tie-gcm
    use fields_module,only: levd0,levd1,lond0,lond1,latd0,latd1
    use mpi_module, only: TIEGCM_WORLD

    ! intern
    use state_module,only: levX1

    implicit none

    ! arguments
    class(sparse_state_interpolator) :: this
    real, dimension(3) :: point ! lon (-180 deg :180 deg), lat(-90 deg : 90 deg), alt (m)
    real, dimension(levd0:levd1,lond0:lond1,latd0:latd1), intent(inout), optional :: &
          zg_mid, zg_mid_nm, zg_int, zg_int_nm
    logical, intent(in), optional :: init_flags(4)
    integer, intent(in), optional :: degree

    ! local

    logical :: init_flags_(4)
    integer :: subdomain_rank

    if( present(init_flags) .eqv. .false.) then
      init_flags_ = .true.
    else
      init_flags_ = init_flags
    end if


    this%on_this_rank = is_within_sub_domain(lon=point(1),lat=point(2))
    subdomain_rank = get_rank_of_subdomain_where_point_is_located(this%on_this_rank)
    if(this%on_this_rank) then
      this%exchange_halo = requires_halo_exchange(lon=point(1),lat=point(2))
    end if
    call MPI_Bcast(this%exchange_halo,1,MPI_LOGICAL,subdomain_rank,TIEGCM_WORLD)

    ! ATTENTION upper TIE-GCM layer contains invalid values. Thus we use 1:levX1

    if( (present(zg_mid).eqv..true.) .and. (init_flags_(1).eqv..true.) )then
      if(cfg_log%verbose_level>0) write(*,*) 'init current midlevel points'
      call this%inter_mid%init(degree=degree,&
                               lon_src=lon_p_halo,&
                               lat_src=lat_p_halo,&
                               lon_dst=point(1),&
                               lat_dst=point(2),&
                               alt_dst=point(3),&
                               subdomain_rank=subdomain_rank)
      if(this%exchange_halo) call fill_halo_and_rim(zg_mid)
      call this%inter_mid%set_geometric_height(zg_mid(1:levX1,:,:))
    end if
    if( (present(zg_mid_nm).eqv..true.) .and. (init_flags_(2).eqv..true.) )then
      if(cfg_log%verbose_level>0) write(*,*) 'init previous midlevel points'
      call this%inter_mid_nm%init(degree=degree,&
                               lon_src=lon_p_halo,&
                               lat_src=lat_p_halo,&
                               lon_dst=point(1),&
                               lat_dst=point(2),&
                               alt_dst=point(3),&
                               subdomain_rank=subdomain_rank)
      if(this%exchange_halo) call fill_halo_and_rim(zg_mid_nm)
      call this%inter_mid_nm%set_geometric_height(zg_mid_nm(1:levX1,:,:))
    end if
    if( (present(zg_int).eqv..true.) .and. (init_flags_(3).eqv..true.) )then
      if(cfg_log%verbose_level>0) write(*,*) 'init current interface points'
      call this%inter_int%init(degree=degree,&
                               lon_src=lon_p_halo,&
                               lat_src=lat_p_halo,&
                               lon_dst=point(1),&
                               lat_dst=point(2),&
                               alt_dst=point(3),&
                               subdomain_rank=subdomain_rank)
      if(this%exchange_halo) call fill_halo_and_rim(zg_int)
      call this%inter_int%set_geometric_height(zg_int(1:levX1,:,:))
    end if
    if( (present(zg_int_nm).eqv..true.) .and. (init_flags_(4).eqv..true.) )then
      if(cfg_log%verbose_level>0) write(*,*) 'init previous interface points'
      call this%inter_int_nm%init(degree=degree,&
                               lon_src=lon_p_halo,&
                               lat_src=lat_p_halo,&
                               lon_dst=point(1),&
                               lat_dst=point(2),&
                               alt_dst=point(3),&
                               subdomain_rank=subdomain_rank)
      if(this%exchange_halo) call fill_halo_and_rim(zg_int_nm)
      call this%inter_int_nm%set_geometric_height(zg_int_nm(1:levX1,:,:))
    end if

    this%prognostic_field_mask = .false.

  end subroutine

  !> Interpolates src_field onto the sparse interpolator's single
  !! destination point, exchanging halo only if required, and marks
  !! prognostic fields.
  subroutine sparse_state_interpolator_interpolate(this, info, src_field, dst, logz)


    ! tie-gcm
    use fields_module,only: levd0,lond0,latd0, f4d

    ! intern
    use array_mapping_module, only: map_to_3d
    use array_print_module, only: printMat
    use character_routines_module, only: string_ends_with
    use quantity_info_module
    use state_module, only: levX1

    implicit none

    ! arguments
    class(sparse_state_interpolator) :: this
    type(quantity_info) :: info
    ! assumed shape in the vertical, so that a quantity without vertical extent
    ! is rejected below instead of being read out of bounds
    real, dimension(levd0:,lond0:,latd0:), intent(inout) :: src_field
    real, intent(out) :: dst
    logical, intent(in), optional :: logz

    ! local

    logical :: logz_

    integer :: idx

    ! LEVEL_NONE is not supported here. Interpolating a quantity without
    ! vertical extent to a single point (e.g. along a satellite track) is
    ! currently not required, thus no horizontal sparse interpolator exists.
    ! This is checked up front, since such a quantity is stored with a single
    ! level only and must not reach the routines below.
    if(info%level==LEVEL_NONE) then
      call shutdown('sparse_state_interpolator_interpolate: quantities without vertical'// &
                    ' extent are not supported. Quantity: '//trim(info%name))
    end if

    if(present(logz)) then
      logz_=logz
    else
      logz_ = .false.
      ! DEN/DEN_NM
      select case (trim(info%name))
        case ("DEN", "DEN_NM")
          logz_ = .true.
      end select
    end if

    ! ATTENTION after interpolation, halo of prognostic fields should be set back to zero,
    ! else the state is slighly changed in the next step. This is done in destoy subroutine
    if(this%exchange_halo)then
      call fill_halo_and_rim(src_field)

      idx = findloc(f4d%short_name,info%name,dim=1)
      if(idx>0) then
        this%prognostic_field_mask(idx)=.true.
      end if
    end if

    select case(info%level)
      case(LEVEL_MID)
        select case(info%step)
          case(STEP_CURRENT)
            call this%inter_mid%interpolate(src_field(levd0:levX1,:,:),dst,logz_)
          case(STEP_PREVIOUS)
            call this%inter_mid_nm%interpolate(src_field(levd0:levX1,:,:),dst,logz_)
        end select
      case(LEVEL_INT)
        select case(info%step)
          case(STEP_CURRENT)
            call this%inter_int%interpolate(src_field(levd0:levX1,:,:),dst,logz_)
          case(STEP_PREVIOUS)
            call this%inter_int_nm%interpolate(src_field(levd0:levX1,:,:),dst,logz_)
        end select
      case default
        ! LEVEL_NONE is not supported here. Interpolating a quantity without
        ! vertical extent to a single point (e.g. along a satellite track) is
        ! currently not required, thus no horizontal sparse interpolator exists.
        call shutdown('sparse_state_interpolator_interpolate: unhandled level of quantity '//trim(info%name))
    end select

!     call printMat(dst_ptr,'interp',order=(/3,2,1/))

  end subroutine

  !> Destroys the sparse interpolator components and zeroes the halo of
  !! any prognostic field that was interpolated.
  subroutine sparse_state_interpolator_destroy(this)

    ! tie-gcm
    use fields_module, only: f4d,itc

    implicit none

    ! arguments
    class(sparse_state_interpolator) :: this

    ! local
    integer :: fidx

    call this%inter_int%destroy
    call this%inter_int_nm%destroy
    call this%inter_mid%destroy
    call this%inter_mid_nm%destroy

    do fidx = 1, size(this%prognostic_field_mask)
      if(this%prognostic_field_mask(fidx).eqv..true.)then
        if(cfg_log%verbose_level>0) write(*,*) f4d(fidx)%short_name, " set halo zero"
        call zero_halo_and_rim(f4d(fidx)%data(:,:,:,itc))
      end if
    end do

  end subroutine

  !> Builds the longitude and latitude B-spline interpolation matrices for
  !! a regular-grid spline interpolator, given source and destination
  !! coordinates.
  subroutine reg_grid_spline_interpolator_init(this, lon_src, lat_src, lon_dst, lat_dst, alt_dst, degree)
    implicit none

    ! arguments
    class(reg_grid_spline_interpolator) :: this

    ! longitudes and latituds of source and destination grid on this rank
    real, dimension(:), intent(in) :: lon_src ! including ghost cells
    real, dimension(:), intent(in) :: lat_src ! including ghost cells
    real, dimension(:), intent(in) :: lon_dst
    real, dimension(:), intent(in) :: lat_dst
    real, dimension(:), intent(in) :: alt_dst

    integer, intent(in), optional :: degree

    if(present(degree))then
      this%degree = degree
    end if


    this%n_lon_src=size(lon_src, dim=1)
    this%n_lat_src=size(lat_src, dim=1)

    allocate( this%alt_dst, source=alt_dst )
    allocate( this%lon_dst, source=lon_dst )
    allocate( this%lat_dst, source=lat_dst )

    call init_horizontal_splines(this%lon_spline, this%lat_spline, &
                                 this%lon_jacobian, this%lat_jacobian, &
                                 lon_src, lat_src, this%degree)

    allocate(this%alt_spline(this%n_lon_src,this%n_lat_src))
    allocate(this%alt_jacobian(this%n_lon_src,this%n_lat_src))

  end subroutine

  !> Builds the longitude/latitude B-spline knots and Jacobians shared by
  !! reg_grid_spline_interpolator and reg_grid_spline_horizontal_interpolator.
  subroutine init_horizontal_splines(lon_spline, lat_spline, lon_jacobian, lat_jacobian, &
                                      lon_src, lat_src, degree)
    implicit none

    ! arguments
    type(bspline), intent(inout) :: lon_spline, lat_spline
    type(sparse_band_matrix), intent(inout) :: lon_jacobian, lat_jacobian
    real, dimension(:), intent(in) :: lon_src ! including ghost cells
    real, dimension(:), intent(in) :: lat_src ! including ghost cells
    integer, intent(in) :: degree

    call lon_spline%compute_knots_for_interpolation(lon_src,degree)
    call lon_spline%compute_jacobian_matrix(lon_src,lon_jacobian)

    call lat_spline%compute_knots_for_interpolation(lat_src,degree)
    call lat_spline%compute_jacobian_matrix(lat_src,lat_jacobian)

  end subroutine

  !> Deallocates all spline and Jacobian data held by the interpolator.
  subroutine reg_grid_spline_interpolator_destroy(this)
    implicit none

    ! arguments
    class(reg_grid_spline_interpolator) :: this

    ! local
    integer :: i,j

    if(allocated(this%alt_dst)) deallocate(this%alt_dst)
    if(allocated(this%lon_dst)) deallocate(this%lon_dst)
    if(allocated(this%lat_dst)) deallocate(this%lat_dst)

    call this%lat_spline%deallocate()
    call this%lon_spline%deallocate()

    call this%lat_jacobian%deallocate()
    call this%lon_jacobian%deallocate()

    if(allocated(this%alt_spline)) then
      do i = 1, size(this%alt_spline,dim=1)
        do j = 1, size(this%alt_spline,dim=2)
            call this%alt_spline(i,j)%deallocate()
        end do
      end do
      deallocate(this%alt_spline)
    end if

    if(allocated(this%alt_jacobian)) then
      do i = 1, size(this%alt_jacobian,dim=1)
        do j = 1, size(this%alt_jacobian,dim=2)
            call this%alt_jacobian(i,j)%deallocate()
        end do
      end do
      deallocate(this%alt_jacobian)
    end if
  end subroutine

  !> Computes the vertical B-spline knots and Jacobians at each source
  !! lon/lat column from the geometric height field Z.
  subroutine reg_grid_spline_interpolator_set_geometric_height(this, Z)

    implicit none

    ! arguments
    class(reg_grid_spline_interpolator) :: this
    real, dimension(:,:,:), intent(in) :: Z

    ! local
    integer :: i,j
    character(len=256) :: errmsg
    real, dimension(size(Z,dim=1)-1,size(Z,dim=2),size(Z,dim=3)) :: dZ
    integer :: loc(3)

    real :: lb, ub

    ! compute_knots_for_interpolation requires Z(:,i,j) to be
    ! strictly increasing with level at every column. An unphysical
    ! model state (e.g. a bad analysis increment) can violate this, and
    ! would otherwise surface only as an opaque out-of-bounds crash deep
    ! inside the spline library (dierckx/fpbspl.f)
    ! TODO check this in gfu lib (compute_knots_for_interpolation)
    ! TODO make this an optional debug option
    if (any(isnan(Z))) then
      write(errmsg,*) 'reg_grid_spline_interpolator_set_geometric_height: ',&
        'geometric height field Z contains NaN'
      call shutdown(trim(errmsg))
    end if

    dZ = Z(2:size(Z,dim=1),:,:) - Z(1:size(Z,dim=1)-1,:,:)
    if (any(dZ <= 0.)) then
      loc = minloc(dZ) ! (level,lon_src,lat_src) index of the worst violation
      write(errmsg,*) 'reg_grid_spline_interpolator_set_geometric_height: ',&
        'Z not strictly increasing at lon_src index=',loc(2),' lat_src index=',loc(3),&
        ' between level',loc(1),'Z=',Z(loc(1),loc(2),loc(3)),&
        ' and level',loc(1)+1,'Z=',Z(loc(1)+1,loc(2),loc(3))
      call shutdown(trim(errmsg))
    end if

    ub = max(MAXVAL(Z(size(Z,dim=1),:,:))/100,this%alt_dst(size(this%alt_dst,dim=1)))
    lb = min(MINVAL(Z(1,:,:))/100,this%alt_dst(1))

    do i = 1,this%n_lon_src
      do j = 1, this%n_lat_src
        call this%alt_spline(i,j)%compute_knots_for_interpolation(Z(:,i,j)/100,this%degree,lb=lb,ub=ub)
        call this%alt_spline(i,j)%compute_jacobian_matrix(Z(:,i,j)/100,this%alt_jacobian(i,j))
      end do
    end do

  end subroutine

  !> Interpolates src onto the regular destination grid via sequential
  !! vertical, longitudinal, and latitudinal B-spline interpolation,
  !! optionally in log space.
  subroutine reg_grid_spline_interpolator_interpolate(this, src, dst, logz)

    ! intern
    use array_print_module, only: printMat
    use netcdf_functionality, only: write_mat_to_netcdf

    implicit none

    ! arguments
    class(reg_grid_spline_interpolator) :: this
    real, dimension(:,:,:), intent(in) :: src  ! (lev,lon,lat)
    real, dimension(:,:,:), intent(out) :: dst ! (lon,lat,alt)
    logical, intent(in) :: logz

    ! local
    integer :: i, j
    character(len=128) :: errmsg

    integer :: n_lon_dst, n_lat_dst, n_alt_dst

    real, dimension(:,:,:), allocatable :: along_alt
    real, dimension(:,:,:), allocatable :: along_lon

    n_lon_dst = size(this%lon_dst,dim=1)
    n_lat_dst = size(this%lat_dst,dim=1)
    n_alt_dst = size(this%alt_dst,dim=1)

    if (size(dst,dim=1) /= n_lon_dst) then
      write(errmsg,*) 'tiegcm_optimized_interpolator: invalid first dimension. ',&
        'nlon of provided dst output: ', size(dst,dim=1), ' expected:', n_lon_dst
      call shutdown(trim(errmsg))
    end if

    if (size(dst,dim=2) /= n_lat_dst) then
      write(errmsg,*) 'tiegcm_optimized_interpolator: invalid second dimension. ',&
        'nlat of provided dst output: ', size(dst,dim=2), ' expected:', n_lat_dst
      call shutdown(trim(errmsg))
    end if

    if (size(dst,dim=3) /= n_alt_dst) then
      write(errmsg,*) 'tiegcm_optimized_interpolator: invalid third dimension. ',&
        'nlalt of provided dst output: ', size(dst,dim=3), ' expected:', n_alt_dst
      call shutdown(trim(errmsg))
    end if

    allocate(along_alt(this%n_lon_src,&
                       this%n_lat_src,&
                       n_alt_dst))

    allocate(along_lon(n_lon_dst,&
                       this%n_lat_src,&
                       n_alt_dst))

    do i = 1, this%n_lon_src
      do j = 1, this%n_lat_src
        if(logz)then
          call this%alt_spline(i,j)%solve(this%alt_jacobian(i,j),log(src(:,i,j)))
          call this%alt_spline(i,j)%eval(this%alt_dst,along_alt(i,j,:))
          along_alt(i,j,:) = exp(along_alt(i,j,:))
        else
          call this%alt_spline(i,j)%solve(this%alt_jacobian(i,j),src(:,i,j))
          call this%alt_spline(i,j)%eval(this%alt_dst,along_alt(i,j,:))
        end if
      end do
    end do

!     call write_mat_to_netcdf(along_alt,"along_alt")

    do i = 1, n_alt_dst
      do j = 1, this%n_lat_src
        call this%lon_spline%solve(this%lon_jacobian,along_alt(:,j,i))
        call this%lon_spline%eval(this%lon_dst,along_lon(:,j,i))
      end do
    end do

!     call write_mat_to_netcdf(along_lon,"along_lon")

    do i = 1, n_alt_dst
      do j = 1, n_lon_dst
        call this%lat_spline%solve(this%lat_jacobian,along_lon(j,:,i))
        call this%lat_spline%eval(this%lat_dst,dst(j,:,i))
      end do
    end do

    deallocate(along_alt,along_lon)

  end subroutine

  !> Initializes the longitude/latitude B-spline stencil for a quantity with
  !! no vertical extent (e.g. VTEC).
  subroutine reg_grid_spline_horizontal_interpolator_init(this, lon_src, lat_src, lon_dst, lat_dst, degree)
    implicit none

    ! arguments
    class(reg_grid_spline_horizontal_interpolator) :: this

    ! longitudes and latitudes of source and destination grid on this rank
    real, dimension(:), intent(in) :: lon_src ! including ghost cells
    real, dimension(:), intent(in) :: lat_src ! including ghost cells
    real, dimension(:), intent(in) :: lon_dst
    real, dimension(:), intent(in) :: lat_dst

    integer, intent(in), optional :: degree

    if(present(degree))then
      this%degree = degree
    end if

    this%n_lon_src=size(lon_src, dim=1)
    this%n_lat_src=size(lat_src, dim=1)

    allocate( this%lon_dst, source=lon_dst )
    allocate( this%lat_dst, source=lat_dst )

    call init_horizontal_splines(this%lon_spline, this%lat_spline, &
                                 this%lon_jacobian, this%lat_jacobian, &
                                 lon_src, lat_src, this%degree)

  end subroutine

  !> Deallocates all spline and Jacobian data held by the interpolator.
  subroutine reg_grid_spline_horizontal_interpolator_destroy(this)
    implicit none

    ! arguments
    class(reg_grid_spline_horizontal_interpolator) :: this

    if(allocated(this%lon_dst)) deallocate(this%lon_dst)
    if(allocated(this%lat_dst)) deallocate(this%lat_dst)

    call this%lat_spline%deallocate()
    call this%lon_spline%deallocate()

    call this%lat_jacobian%deallocate()
    call this%lon_jacobian%deallocate()

  end subroutine

  !> Interpolates a quantity with no vertical extent
  !! horizontally onto the regular destination grid,
  !! via sequential longitude then latitude B-spline interpolation.
  !! Both src and dst carry a size-1 vertical axis for compatibility
  !! with interfaces especting 3D fields
  subroutine reg_grid_spline_horizontal_interpolator_interpolate(this, src, dst)

    implicit none

    ! arguments
    class(reg_grid_spline_horizontal_interpolator) :: this
    real, dimension(:,:,:), intent(in) :: src  ! (1,lon,lat)
    real, dimension(:,:,:), intent(out) :: dst ! (lon,lat,1)

    ! local
    integer :: j
    character(len=128) :: errmsg

    integer :: n_lon_dst, n_lat_dst

    real, dimension(:,:), allocatable :: along_lon

    n_lon_dst = size(this%lon_dst,dim=1)
    n_lat_dst = size(this%lat_dst,dim=1)

    if (size(src,dim=1) /= 1) then
      write(errmsg,*) 'tiegcm_optimized_interpolator: reg_grid_spline_horizontal_interpolator ',&
        'expects a source field with no vertical extent (first dimension size 1), got: ', size(src,dim=1)
      call shutdown(trim(errmsg))
    end if

    if (size(dst,dim=1) /= n_lon_dst) then
      write(errmsg,*) 'tiegcm_optimized_interpolator: invalid first dimension. ',&
        'nlon of provided dst output: ', size(dst,dim=1), ' expected:', n_lon_dst
      call shutdown(trim(errmsg))
    end if

    if (size(dst,dim=2) /= n_lat_dst) then
      write(errmsg,*) 'tiegcm_optimized_interpolator: invalid second dimension. ',&
        'nlat of provided dst output: ', size(dst,dim=2), ' expected:', n_lat_dst
      call shutdown(trim(errmsg))
    end if

    if (size(dst,dim=3) /= 1) then
      write(errmsg,*) 'tiegcm_optimized_interpolator: invalid third dimension. ',&
        'nalt of provided dst output: ', size(dst,dim=3), ' expected: 1'
      call shutdown(trim(errmsg))
    end if

    allocate(along_lon(n_lon_dst, this%n_lat_src))

    do j = 1, this%n_lat_src
      call this%lon_spline%solve(this%lon_jacobian,src(1,:,j))
      call this%lon_spline%eval(this%lon_dst,along_lon(:,j))
    end do

    do j = 1, n_lon_dst
      call this%lat_spline%solve(this%lat_jacobian,along_lon(j,:))
      call this%lat_spline%eval(this%lat_dst,dst(j,:,1))
    end do

    deallocate(along_lon)

  end subroutine

  !> Initializes the longitude/latitude B-spline stencil for a single
  !! destination point, on the point's owning subdomain rank only.
  subroutine sparse_interpolator_init(this,lon_src,lat_src,lon_dst,lat_dst,alt_dst,degree,subdomain_rank)

    use mpi_module, only: mytid

    implicit none

    ! arguments
    class(sparse_interpolator) :: this
    ! longitudes and latituds of source grid on this rank
    real, dimension(:), intent(in) :: lon_src ! including ghost cells
    real, dimension(:), intent(in) :: lat_src ! including ghost cells
    real, intent(in) :: lon_dst
    real, intent(in) :: lat_dst
    real, intent(in) :: alt_dst
    integer, intent(in), optional :: degree
    integer, intent(in), optional :: subdomain_rank

    if(present(degree))then
      if (degree>3) then
        write(*,*) 'maximal spline degree is 3 (2 ghost cells at borders)'
        this%degree = 3
      else
        this%degree= degree
      end if
    end if

    this%lon_dst = lon_dst
    this%lat_dst = lat_dst
    this%alt_dst = alt_dst

    if(present(subdomain_rank))then
      this%on_this_rank = (mytid==subdomain_rank)
      this%subdomain_rank = subdomain_rank
    else
      this%on_this_rank = is_within_sub_domain(lon_dst,lat_dst)
      this%subdomain_rank = get_rank_of_subdomain_where_point_is_located(this%on_this_rank)
    end if


    if( this%on_this_rank ) then

      ! TODO Works also for for splines with degree smaller than 3, but perfoms unessesary computations
      call get_range_for_non_zero_cubic_bsplines(lon_src, lon_dst, this%lon_first, this%lon_last)
      call get_range_for_non_zero_cubic_bsplines(lat_src, lat_dst, this%lat_first, this%lat_last)

      allocate(this%alt_spline(4,4))
      allocate(this%alt_jacobian(4,4))

      call this%lon_spline%compute_knots_for_interpolation(lon_src(this%lon_first:this%lon_last),&
                                                          this%degree)
      call this%lon_spline%compute_jacobian_matrix(lon_src(this%lon_first:this%lon_last),&
                                                  this%lon_jacobian)

      call this%lat_spline%compute_knots_for_interpolation(lat_src(this%lat_first:this%lat_last),&
                                                          this%degree)
      call this%lat_spline%compute_jacobian_matrix(lat_src(this%lat_first:this%lat_last),&
                                                  this%lat_jacobian)
     end if
  end subroutine

  !> Deallocates the spline and Jacobian data of a sparse_interpolator.
  subroutine sparse_interpolator_destroy(this)
    implicit none

    ! arguments
    class(sparse_interpolator) :: this

    ! local
    integer :: i,j

    call this%lat_spline%deallocate()
    call this%lon_spline%deallocate()

    call this%lat_jacobian%deallocate()
    call this%lon_jacobian%deallocate()

    if(allocated(this%alt_spline)) then
      do i = 1, size(this%alt_spline,dim=1)
        do j = 1, size(this%alt_spline,dim=2)
            call this%alt_spline(i,j)%deallocate()
        end do
      end do
      deallocate(this%alt_spline)
    end if

    if(allocated(this%alt_jacobian)) then
      do i = 1, size(this%alt_jacobian,dim=1)
        do j = 1, size(this%alt_jacobian,dim=2)
            call this%alt_jacobian(i,j)%deallocate()
        end do
      end do
      deallocate(this%alt_jacobian)
    end if
  end subroutine

  !> Computes the vertical B-spline knots and Jacobians for the local
  !! stencil, on the point's owning rank only.
  subroutine sparse_interpolator_set_geometric_height(this, z)

    implicit none
    ! arguments
    class(sparse_interpolator) :: this
    real, dimension(:,:,:), intent(in) :: Z ! including ghost cells

    ! local
    integer :: m
    real, dimension(:,:,:), allocatable :: Zsub
    real :: lb, ub
    integer :: i,j

    if( this%on_this_rank ) then
      m = size(this%alt_spline,dim=1)

      allocate(Zsub,source=Z(:,this%lon_first:this%lon_last,this%lat_first:this%lat_last))

      ub = max(MAXVAL(Zsub(size(Zsub,dim=1),:,:))/100,this%alt_dst)
      lb = min(MINVAL(Zsub(1,:,:))/100,this%alt_dst)

      do i = 1, m
        do j = 1, m
          call this%alt_spline(i,j)%compute_knots_for_interpolation(Zsub(:,i,j)/100,this%degree,lb=lb,ub=ub)
          call this%alt_spline(i,j)%compute_jacobian_matrix(Zsub(:,i,j)/100,this%alt_jacobian(i,j))
        end do
      end do

      deallocate(Zsub)
    end if

  end subroutine

  !> Interpolates src to the interpolator's single destination point on
  !! the owning rank, then broadcasts the result to all ranks.
  subroutine sparse_interpolator_interpolate(this, src, dst, logz)

    ! intern
    use array_print_module, only: printMat
    use netcdf_functionality, only: write_mat_to_netcdf

    implicit none

    ! arguments
    class(sparse_interpolator) :: this
    real, dimension(:,:,:), intent(in) :: src  ! (lev,lon,lat) including ghost cells
    real, intent(out) :: dst
    logical, intent(in) :: logz

    ! local
    integer :: i, j

    integer  :: m

    real, dimension(:,:), allocatable :: along_alt
    real, dimension(:), allocatable :: along_lon

    real, dimension(:,:,:), allocatable :: src_sub

    if( this%on_this_rank ) then
      m =4

      allocate(src_sub,source=src(:,this%lon_first:this%lon_last,this%lat_first:this%lat_last))

      allocate(along_alt(m,m))

      allocate(along_lon(m))

      do i = 1, m
        do j = 1, m
          if(logz)then
            call this%alt_spline(i,j)%solve(this%alt_jacobian(i,j),log(src_sub(:,i,j)))
            call this%alt_spline(i,j)%eval(this%alt_dst,along_alt(i,j))
            along_alt(i,j) = exp(along_alt(i,j))
          else
            call this%alt_spline(i,j)%solve(this%alt_jacobian(i,j),src_sub(:,i,j))
            call this%alt_spline(i,j)%eval(this%alt_dst,along_alt(i,j))
          end if
        end do
      end do

  !     call write_mat_to_netcdf(along_alt,"along_alt")

      do j = 1, m
        call this%lon_spline%solve(this%lon_jacobian,along_alt(:,j))
        call this%lon_spline%eval(this%lon_dst,along_lon(j))
      end do


  !     call write_mat_to_netcdf(along_lon,"along_lon")

      call this%lat_spline%solve(this%lat_jacobian,along_lon(:))
      call this%lat_spline%eval(this%lat_dst,dst)

      deallocate(along_alt,along_lon,src_sub)
    end if

    call broadcast_interpolated_value_to_all_subdomians(this, dst)

  end subroutine

  !> Determines and broadcasts which single MPI rank owns a point, from
  !! each rank's on_this_rank flag.
  function get_rank_of_subdomain_where_point_is_located(on_this_rank) result(rank)

    ! extern
    use mpi_f08

    ! intern
    use mod_parallel_pdaf, only: local_ntask
    use mpi_module, only: TIEGCM_WORLD, mytid, handle_mpi_err

    implicit none

    ! arguments
    logical, intent(in) :: on_this_rank

    ! result
    integer :: rank

    ! local
    logical, dimension(local_ntask) :: rank_mask
    integer :: ier


    call mpi_gather(on_this_rank,1,MPI_LOGICAL,&
                    rank_mask,1,MPI_LOGICAL,&
                    0, TIEGCM_WORLD, ier)
    if (ier /= 0) call handle_mpi_err(ier,'get_rank_of_subdomain_where_point_is_located')

    if(mytid==0)then
      if(count(rank_mask)==1) then
        rank = findloc(rank_mask,.true.,dim=1)-1
      elseif(count(rank_mask)==0) then
        call shutdown('on_this_rank was not determined')
      else
        call shutdown('on_this_rank is ambiguous')
      end if
    end if

    call mpi_bcast(rank,1,MPI_INT,0,TIEGCM_WORLD,ier)
    if (ier /= 0) call handle_mpi_err(ier,'get_rank_of_subdomain_where_point_is_located')

  end function

  !> Broadcasts an interpolated value from the interpolator's owning
  !! subdomain rank to all ranks.
  subroutine broadcast_interpolated_value_to_all_subdomians(this, val)

    ! extern
    use mpi_f08

    ! intern
    use mpi_module, only: TIEGCM_WORLD, handle_mpi_err

    implicit none

    ! arguments
    class(sparse_interpolator) :: this
    real, intent(inout) :: val

    ! local
    integer :: ier
    call mpi_bcast(val,1,MPI_REAL8,this%subdomain_rank,TIEGCM_WORLD,ier)
    if (ier /= 0) call handle_mpi_err(ier,'broadcast_interpolated_value_to_all_subdomians')

  end subroutine

  !> Debug/test routine: interpolates a synthetic field with
  !! reg_grid_spline_interpolator and writes the result and ground truth
  !! to NetCDF, then aborts.
  subroutine test_reg_grid_interpolator

    ! exterm
    use mpi_f08

    ! tie-gcm
    use fields_module, only: zg,levd0,levd1,lond0,lond1,latd0,latd1
    use mpi_module, only: mytid

    ! intern
    use netcdf_functionality, only: write_mat_to_netcdf
    use test_values, only: synthetic_val, synthetic_val_lin_fun

    implicit none

    type(reg_grid_spline_interpolator) :: inter

    real, dimension(:,:,:), allocatable, target :: src
    real, dimension(:,:,:), allocatable, target :: dst
    real, dimension(:,:,:), allocatable :: truth

    real, dimension(:), allocatable :: lon_dst_all, lat_dst_all, alt_dst_all

    integer :: i,j,k

    integer :: count_lon, offset_lon, count_lat, offset_lat


    call fill_halo_and_rim(zg)

!      call write_mat_to_netcdf(zg,"zg")

    allocate(lon_dst_all(72))
    allocate(lat_dst_all(36))
    allocate(alt_dst_all(5))
    do i=1,size(lon_dst_all,dim=1)
      lon_dst_all(i) = -180+(360./size(lon_dst_all,dim=1))*(i-1)
    end do
    do i=1,size(lat_dst_all,dim=1)
      lat_dst_all(i) = -80+(160./size(lat_dst_all,dim=1))*(i-1)
    end do
    do i=1,size(alt_dst_all,dim=1)
       alt_dst_all(i) = (100+(200/size(alt_dst_all,dim=1))*(i-1))*1000 ! m
    end do

    call distriubute_grid(lon_dst_all, lat_dst_all, offset_lon, count_lon, offset_lat, count_lat)

    allocate(src(levd0:levd1,lond0:lond1,latd0:latd1))
    allocate(dst(count_lon,count_lat,size(alt_dst_all,dim=1)))
    allocate(truth(count_lon,count_lat,size(alt_dst_all,dim=1)))

    ! is not periodic in longitude ! Expect errors when interpolating on periodic boundary
    do i = levd0,levd1
      do j = lond0,lond1
        do k = latd0,latd1
          src(i,j,k) = synthetic_val( lon_p_halo(j), lat_p_halo(k), zg(i,j,k)/100 ) + 1e6
        end do
      end do
    end do

!     call write_mat_to_netcdf(src,"src")

    do i = 1,size(alt_dst_all,dim=1)
      do j = 1,count_lon
        do k = 1,count_lat
          truth(j,k,i) = synthetic_val( lon_dst_all(offset_lon-1+j), lat_dst_all(offset_lat-1+k), alt_dst_all(i) ) + 1e6
        end do
      end do
    end do

    write(*,*) 'rank ', mytid, ' lons:', lon_dst_all(offset_lon:offset_lon+count_lon-1)
    write(*,*) 'rank ', mytid, ' lats:', lat_dst_all(offset_lat:offset_lat+count_lat-1)

    call inter%init(lon_p_halo,lat_p_halo, lon_dst_all(offset_lon:offset_lon+count_lon-1),&
                                  lat_dst_all(offset_lat:offset_lat+count_lat-1),&
                                  alt_dst_all, 3)

    call inter%set_geometric_height(zg)

    call inter%interpolate(src,dst,.false.)

    call write_mat_to_netcdf(truth,"truth")
    call write_mat_to_netcdf(dst,"dst")

    call MPI_barrier(MPI_COMM_WORLD)
    call shutdown('end of test')

    deallocate(src)
    deallocate(dst)
    deallocate(truth)
    deallocate(lon_dst_all)
    deallocate(lat_dst_all)
    deallocate(alt_dst_all)
    call inter%destroy()

  end subroutine

  !> Debug/test routine: interpolates a synthetic field with no vertical
  !! extent with reg_grid_spline_horizontal_interpolator and writes the
  !! result and ground truth to NetCDF, then aborts.
  subroutine test_reg_grid_spline_horizontal_interpolator

    ! exterm
    use mpi_f08

    ! tie-gcm
    use fields_module, only: lond0,lond1,latd0,latd1
    use mpi_module, only: mytid

    ! intern
    use netcdf_functionality, only: write_mat_to_netcdf
    use test_values, only: synthetic_val,synthetic_val_lin_fun

    implicit none

    type(reg_grid_spline_horizontal_interpolator) :: inter

    real, dimension(:,:,:), allocatable, target :: src
    real, dimension(:,:,:), allocatable, target :: dst
    real, dimension(:,:,:), allocatable :: truth

    real, dimension(:), allocatable :: lon_dst_all, lat_dst_all

    ! The field has no vertical extent. This height only fixes the horizontal
    ! pattern of the synthetic field (compare the single shell height of a
    ! VTEC map).
    real, parameter :: shell_height = 450e+3 ! m

    integer :: i,j,k

    integer :: count_lon, offset_lon, count_lat, offset_lat

    allocate(lon_dst_all(72))
    allocate(lat_dst_all(36))
    do i=1,size(lon_dst_all,dim=1)
      lon_dst_all(i) = -180+(360./size(lon_dst_all,dim=1))*(i-1)
    end do
    do i=1,size(lat_dst_all,dim=1)
      lat_dst_all(i) = -80+(160./size(lat_dst_all,dim=1))*(i-1)
    end do

    call distriubute_grid(lon_dst_all, lat_dst_all, offset_lon, count_lon, offset_lat, count_lat)

    allocate(src(1,lond0:lond1,latd0:latd1))
    allocate(dst(count_lon,count_lat,1))
    allocate(truth(count_lon,count_lat,1))

    ! is not periodic in longitude ! Expect errors when interpolating on periodic boundary
    do j = lond0,lond1
      do k = latd0,latd1
        src(1,j,k) = synthetic_val_lin_fun( lon_p_halo(j), lat_p_halo(k), shell_height ) + 1e6
      end do
    end do

!     call write_mat_to_netcdf(src,"src")

    do j = 1,count_lon
      do k = 1,count_lat
        truth(j,k,1) = synthetic_val_lin_fun( lon_dst_all(offset_lon-1+j), lat_dst_all(offset_lat-1+k), shell_height ) + 1e6
      end do
    end do

    write(*,*) 'rank ', mytid, ' lons:', lon_dst_all(offset_lon:offset_lon+count_lon-1)
    write(*,*) 'rank ', mytid, ' lats:', lat_dst_all(offset_lat:offset_lat+count_lat-1)

    call inter%init(lon_p_halo,lat_p_halo, lon_dst_all(offset_lon:offset_lon+count_lon-1),&
                                  lat_dst_all(offset_lat:offset_lat+count_lat-1),&
                                  3)

    call inter%interpolate(src,dst)

    call write_mat_to_netcdf(truth,"truth")
    call write_mat_to_netcdf(dst,"dst")

    call MPI_barrier(MPI_COMM_WORLD)
    call shutdown('end of test')

    deallocate(src)
    deallocate(dst)
    deallocate(truth)
    deallocate(lon_dst_all)
    deallocate(lat_dst_all)
    call inter%destroy()

  end subroutine

  !> Returns the index of the nearest sorted-array element less than or
  !! equal to val (index 1 if val is below the array, n-1 if at or above
  !! the top).
  function find_index_0(array, val) result(idx0)

    ! intern
    use search_module, only: lower_bound

    implicit none

    ! arguments
    real, dimension(:), intent(in) :: array ! sorted coordinates
    real, intent(in) :: val ! queried value

    ! local
    integer :: idx0
    integer :: n

    n = ubound(array,dim=1)

    if(val >= array(n)) then
      idx0 = n-1

    else if (val <= array(1)) then
      idx0 = 1

    else
      idx0 = lower_bound(array,val)
      if(array(idx0)/=val)then
        idx0=idx0-1
      end if
    end if

  end function find_index_0

  !> Determines the contiguous slice of a global lon/lat grid that falls
  !! within this rank's subdomain (including ghost cells).
  subroutine distriubute_grid(grid_lons,grid_lats, offset_lon, count_lon, offset_lat, count_lat)

    ! tie-gcm
    use mpi_module, only: mytid

    ! intern
    use search_module, only: lower_bound, upper_bound

    use state_module, only: lonX0, lonX1, latX0, latX1

    implicit none

    ! arguments
    real, dimension(:), intent(in) :: grid_lons ! gathered longitudes of grid
    real, dimension(:), intent(in) :: grid_lats ! gathered latitudes of grid

    integer, intent(out) :: offset_lon
    integer, intent(out) :: count_lon

    integer, intent(out) :: offset_lat
    integer, intent(out) :: count_lat

    ! local
    integer :: idx0,idx1
    integer :: grid_nlon, grid_nlat

    integer :: last_lon, last_lat
    integer :: first_lon, first_lat

    last_lon = ubound(lon_p_halo,dim=1)
    last_lat = ubound(lat_p_halo,dim=1)
    first_lon = lbound(lon_p_halo,dim=1)
    first_lat = lbound(lat_p_halo,dim=1)

    count_lat = 0
    count_lon = 0
    offset_lon = 0
    offset_lat = 0

    grid_nlon = size(grid_lons,dim=1)
    grid_nlat = size(grid_lats,dim=1)

    if(lon_p_halo(first_lon)>grid_lons(grid_nlon)) then
      write(*,*) "this task does not contain any elementes located on observation grid"
      return
    end if
    if(lon_p_halo(last_lon)<grid_lons(1)) then
      write(*,*) "this task does not contain any elementes located on observation grid"
      return
    end if

    if(lat_p_halo(first_lat)>grid_lats(grid_nlat)) then
      write(*,*) "this task does not contain any elementes located on observation grid"
      return
    end if
    if(lat_p_halo(last_lat)<grid_lats(1)) then
      write(*,*) "this task does not contain any elementes located on observation grid"
      return
    end if

    ! first element in grid_lons not less than first non ghost/periodic cell
    idx0 = lower_bound(grid_lons,lon_p_halo(lonX0))

    idx1 = lower_bound(grid_lons,lon_p_halo(lonX1+1))

    if(idx1>grid_nlon)then
      ! last value of lon_p_halo was not found in grid_lons
      idx1 = grid_nlon
    else
      idx1=idx1-1
    end if

    offset_lon = idx0
    count_lon = (idx1-idx0)+1

    write(*,'(a,i4,a,f7.1,a,f7.1)') 'longitudes available at rank ', mytid, ' via ghost cells ', &
      lon_p_halo(first_lon),          " " , lon_p_halo(last_lon)
    write(*,'(a,i4,a,f7.1,a,f7.1,a,i4,a)') 'longitudes processed at rank ',mytid,'                 ', &
      grid_lons(idx0), " " , grid_lons(idx1), "(", count_lon, ")"


    ! first element in grid_lats greater than  non ghost/periodic cell
    idx0 = lower_bound(grid_lats,lat_p_halo(latX0))

    idx1 = lower_bound(grid_lats,lat_p_halo(latX1+1))

    ! The subdomains have to cover the whole observation grid, otherwise the
    ! gathered result has gaps. An observation grid may however extend beyond
    ! the model grid, which reaches from -87.5 to 87.5 degrees at the default
    ! resolution. The two clamps below assign such observations to the
    ! southernmost respectively northernmost rank, which are the only ones
    ! that can take them. Those observations are still located within the
    ! ghost cells of that rank, hence they are interpolated and not
    ! extrapolated.
    ! ATTENTION the ghost cells beyond the pole do not continue the field.
    ! They hold the mean over all longitudes surrounding the pole (see
    ! mp_poles, respectively mp_poles_2d for quantities without vertical
    ! extent). Observations beyond the pole are therefore drawn towards that
    ! zonal mean.

    ! southern clamp: observations south of the first model latitude
    if(at_south_end)then
      idx0 = 1
    end if

    ! northern clamp: observations north of the last model latitude. Here the
    ! search for the upper index runs beyond the observation grid instead,
    ! since lat_p_halo(latX1+1) exceeds all of its latitudes.
    if(idx1>grid_nlat)then
      idx1 = grid_nlat
    else
      ! idx1 is the first latitude of the next rank
      idx1=idx1-1
    end if

    offset_lat = idx0
    count_lat = (idx1-idx0)+1

    write(*,'(a,i4,a,f7.1,a,f7.1)') ' latitudes available at rank ', mytid, ' via ghost cells ', &
      lat_p_halo(first_lat),          " " , lat_p_halo(last_lat)
    write(*,'(a,i4,a,f7.1,a,f7.1,a,i4,a)') ' latitudes processed at rank ',mytid,'                 ', &
      grid_lats(idx0), " " , grid_lats(idx1), "(", count_lat, ")"

  end subroutine

  !> Zeroes the two outermost halo/rim cells of field in longitude and
  !! latitude.
  subroutine zero_halo_and_rim(field)

    ! tie-gcm
    use fields_module, only: levd0,levd1,lond0,lond1,latd0,latd1

    ! intern
    use array_print_module, only: printMat

    real, intent(inout) :: field (levd0:levd1,lond0:lond1,latd0:latd1)

    field(:,lond0:lond0+1,:) = 0
    field(:,lond1-1:lond1,:) = 0
    field(:,:,latd0:latd0+1) = 0
    field(:,:,latd1-1:latd1) = 0

!     call printMat(field,"after halo zero")

  end subroutine

  !> Returns the 4-point index range of source coordinates with non-zero
  !! cubic B-spline support around x_dst.
  subroutine get_range_for_non_zero_cubic_bsplines(x_src,x_dst,first,last)
    implicit none

    ! arguments
    real, dimension(:), intent(in) :: x_src
    real, intent(in) :: x_dst
    integer, intent(out) :: first
    integer, intent(out) :: last

    ! local
    integer :: left_knot

    left_knot = find_index_0(x_src,x_dst)
    first = left_knot-1
    last = left_knot+2

  end subroutine

  !> Performs MPI halo exchange and pole averaging on field, then fills
  !! the periodic rim cells at the east/west subdomain ends.
  subroutine fill_halo_and_rim(field)

    ! tie-gcm
    use fields_module, only: levd0,levd1,lond0,lond1,latd0,latd1
    use mpi_module, only: mp_halo_exchanage, mp_poles

    ! intern
    use array_print_module, only: printMat
!     use state_module, only: lonX0, lonX1, latX0, latX1

    implicit none

    ! arguments
    real, intent(inout) :: field (levd0:levd1,lond0:lond1,latd0:latd1)

    ! local
    integer :: nlon, nlat, nalt

    nlon = (lond1-lond0)+1
    nlat = (latd1-latd0)+1
    nalt = (levd1-levd0)+1

    call mp_halo_exchanage(field)

    ! compute halo at poles from mean of all longitudes around pole
    call mp_poles(field)

    ! compute halo at poles by extending each longitude
!       if(at_south_end) then
!         field(:,:,latd0) = field(:,:,latd0+2)
!         field(:,:,latd0+1) = field(:,:,latd0+2)
!       end if
!
!       if(at_north_end) then
!         field(:,:,latd1) = field(:,:,latd1-2)
!         field(:,:,latd1-1) = field(:,:,latd1-2)
!       end if

    ! periodic cells (actually halo makes them obsolete and they are annoying)
    if(at_west_end) then
      field(:,lond0,:) = field(:,lond0+2,:)
      field(:,lond0+1,:) = field(:,lond0+2,:)
    end if

    if(at_east_end) then
      field(:,lond1,:) = field(:,lond1-2,:)
      field(:,lond1-1,:) = field(:,lond1-2,:)
    end if

!     call printMat(field,'field after halo and pole fill')

  end subroutine

  !> Fills the halo/rim of a 2D (lon,lat) field with no vertical extent
  !! (e.g. VTEC). Same as fill_halo_and_rim, but for a field without a
  !! vertical dimension.
  subroutine fill_halo_and_rim_2d(field)

    ! tie-gcm
    use fields_module, only: lond0,lond1,latd0,latd1
    use mpi_module, only: mp_halo_exchanage_2d, mp_poles_2d

    implicit none

    ! arguments
    real, intent(inout) :: field (lond0:lond1,latd0:latd1)

    call mp_halo_exchanage_2d(field)

    ! compute halo at poles from mean of all longitudes around pole
    call mp_poles_2d(field)

    ! periodic cells (actually halo makes them obsolete and they are annoying)
    if(at_west_end) then
      field(lond0,:) = field(lond0+2,:)
      field(lond0+1,:) = field(lond0+2,:)
    end if

    if(at_east_end) then
      field(lond1,:) = field(lond1-2,:)
      field(lond1-1,:) = field(lond1-2,:)
    end if

  end subroutine

  !> Elementwise test of whether each (lon,lat) pair lies within this
  !! rank's subdomain bounding box.
  function is_within_sub_domain_array(lon,lat) result(is_within)

    use mpi_module, only: mytid

    implicit none

    ! arguments
    real, dimension(:), intent(in) :: lon
    real, dimension(size(lon,dim=1)), intent(in) :: lat

    ! result
    logical, dimension(size(lon,dim=1)) :: is_within

    where((lb_lon(mytid) < lon) .and. (lon <= ub_lon(mytid)) .and. &
          (lb_lat(mytid) < lat) .and. (lat <= ub_lat(mytid)))
      is_within = .true.
    else where
      is_within = .false.
    end where

  end function

  !> Tests whether a single (lon,lat) point lies within this rank's
  !! subdomain bounding box.
  function is_within_sub_domain_scalar(lon,lat) result(is_within)

    use mpi_module, only: mytid

    implicit none

    ! arguments
    real, intent(in) :: lon
    real, intent(in) :: lat

    ! result
    logical :: is_within

    is_within = ( (lb_lon(mytid) < lon) .and. (lon <= ub_lon(mytid)) .and. &
                  (lb_lat(mytid) < lat) .and. (lat <= ub_lat(mytid)) )

!     write(*,*) 'is_within:', is_within
!     write(*,*) 'lon: ', lb_lon(mytid), lon, ub_lon(mytid)
!     write(*,*) 'lat: ', lb_lat(mytid), lat, ub_lat(mytid)

  end function

  !> Debug check: gathers on_rank flags on rank 0 and reports whether
  !! exactly one subdomain claims ownership of a point.
  subroutine check_unambigous_sub_domain_assignment(on_rank)

    ! extern
    use mpi_f08

    ! tie-gcm
    use mpi_module, only: mytid, TIEGCM_WORLD

    ! intern
    use mod_parallel_pdaf, only: local_ntask

    implicit none

    ! arguments
    logical, intent(in) :: on_rank

    ! local
    logical, dimension(local_ntask) :: obs_on_rank

    call mpi_gather(on_rank,1,MPI_LOGICAL,&
                    obs_on_rank,1,MPI_LOGICAL,&
                    0,TIEGCM_WORLD)

    if(mytid==0)then
      write(*,*) 'obs domain assignment: ', obs_on_rank
      if( count(obs_on_rank) == 1 ) then
        write(*,*) 'obs was assigned successfully to a subdomian'
      else if( count(obs_on_rank) > 1 ) then
        write(*,*) 'something went wrong, observation was not assigned without ambiguity to a sub domain'
      else
        write(*,*) 'something went wrong, observation was not assigned to any sub domain'
      end if
    end if

  end subroutine

  !> Tests whether a point is close enough to this rank's subdomain
  !! boundary that halo exchange is required before interpolating.
  function requires_halo_exchange(lon,lat) result(do_halo)

    ! intern
    use state_module,only: latX0, latX1, lonX0, lonX1

    implicit none

    ! arguments
    real, intent(in) :: lon
    real, intent(in) :: lat

    ! result
    logical :: do_halo

    ! the fields contain values between [lon_p_halo(lonX0),lon_p_halo(lonX1)] and
    ! [lat_p_halo(latX0),lat_p_halo(latX1)] without halo exchange.
    ! As soon as we are two cells away from that boundary we do not need to perform the exchange

    do_halo = ( (lon <= lon_p_halo(lonX0+2)) .or. (lon >= lon_p_halo(lonX1-2)) .or. &
                (lat <= lat_p_halo(latX0+2)) .or. (lat >= lat_p_halo(latX1-2)) )

!     write(*,*) lon, '[',lon_p_halo(lonX0+2),',', lon_p_halo(lonX1-2),']'
!     write(*,*) lat, '[',lat_p_halo(latX0+2),',', lat_p_halo(latX1-2),']'
!     write(*,*) 'halo exchange necessary', do_halo

  end function

END MODULE tiegcm_optimized_interpolator
