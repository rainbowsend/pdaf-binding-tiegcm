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
! Maps TIE-GCM prognostic fields and calibration parameters into/out of the PDAF state vector, including domain-decomposition bookkeeping.

! This modules stores where which parameters are stored in the state vector.
! The data of the state vectot is not stored in this module. Instead it is an
! argument in the corresponding subroutines
!
! TIEGCM has several fields, e.g, 'TN', 'O1', or 'HE'. This module
! maps the content of the selected fields to the state vector,
! considering the domain decomposition of TIEGCM:
!    The fields (girds) are divided in lat lon parts
!
!   dimensions in tiegcm : lev x lon  x lat [ x time ]
!
!   example of decomposition of a single field using 4 tasks
!           -------------
!          |   3  |   4  |
!      lat  -------------
!          |   1  |   2  |
!           -------------
!              lon
! 
! We need to know where which subdomain of which field is saved
! You can order either by field (F) or rank (R) (see Example below)
!
! order by field        ordered by rank
!
!      | TN  1           | 1 TN
!      | TN  2           | 1 O1
!      | TN  3           | 2 TN
!      | TN  4           | 2 O1
!      | O1  1           | 3 TN
!      | O1  2           | 3 O1
!      | O1  3           | 4 TN
!      | O1  4           | 4 O1
!
! Both methods are useful. F-order is useful for reading hist fields
! and ensemble generation (function 'PDAF_eofcovar'). R-order is useful
! for distributing the state vector with MPI.
!
! be aware! R related indices start at 0, but F-related indices start at 1
!
!
! Armin Corbin
! University of Bonn, APMG
! 18.12.2019
!
!
!  ATTENTION Do not use observation_module within this file! This may lead to circular dependencies
!
module state_module

  ! intern
  use mapping_module, only: construct_mapping, print_mapping, mapping, deconstruct_mapping

  IMPLICIT NONE

  ! Augumented state vector at rank 0
  !
  ! z = |x|  'state'
  !     |p|  'dynamics'
  !
  !
  ! The augumented state vector includes n fields from f4d and m calibration parameters
  !
  ! id in map   content
  !
  !     1       fd_idx(1)     |x_1| <-- idx_f3d_0
  !     2       fd_idx(2)     |x_2|
  !                           |...|
  !     n       fd_idx(n)     |x_n| <-- idx_f3d_1
  !     n+1     cal_idx(1)    |c_1| <-- idx_cal_0
  !     n+2     cal_idx(2)    |c_2|
  !                           |...|
  !     n+m     cal_idx(m)    |c_m| <-- idx_cal_1
  !

  type :: state_vector_mapping
    type(mapping) :: map ! mapping of TIE-GCM fields to the state vector
    integer, dimension(:), allocatable :: fd_idx    ! corresponding idx in f4d from fields.F
    integer, dimension(:), allocatable :: cal_idx    ! corresponding idx in model_parameters from model_parameters_handling.F90
    integer :: idx_f3d_0, idx_f3d_1 ! first type saved in state vector are 3d TIE-GCM fields
    integer :: idx_cal_0, idx_cal_1 ! second type saved in state vector are calibratio parameters
    character(len=16), dimension(:), allocatable :: tgcm_field_names
    character(len=16), dimension(:), allocatable :: calibration_names
    contains
    procedure, pass(this), public :: init => state_vector_mapping_init
    procedure, pass(this), public :: destroy => state_vector_mapping_destroy
    procedure, pass(this), public :: fill_state_p => state_vector_mapping_fill_state_p
    procedure, pass(this), public :: fill_state_p_calibration => state_vector_mapping_fill_state_p_calibration
    procedure, pass(this), public :: distribute_calibration => state_vector_mapping_distribute_calibration
    procedure, pass(this), public :: field_3d => state_vector_mapping_get_field_3d
    procedure, pass(this), public :: field => state_vector_mapping_get_field
  end type

  type(state_vector_mapping), protected :: state_vector

  ! this values are valid for any instance of a state vector ---

  integer, protected :: levx0, levx1, nlevx   ! (pressure) levels used in state vector

  ! holds same information as idx_intern%(mytid)%
  integer, protected :: latX0, latX1, nlonX
  integer, protected :: lonX0, lonX1, nlatX

  type task_specific_dim
    integer :: nlats     ! number of latitudes  calculated by this task
    integer :: nlons     ! number of longitudes calculated by this task
    integer :: lat0,lat1 ! first and last latitude  indices
    integer :: lon0,lon1 ! first and last longitude indices
  end type task_specific_dim

  type(task_specific_dim), allocatable, protected :: idx_intern(:), idx_nc(:)

  real, allocatable, dimension(:), protected :: lon_p, lat_p ! longitudes and latitudes on current subdomain

  contains

  !> Computes this process's local level/longitude/latitude index bounds for the
  !! state vector from TIE-GCM's task decomposition.
  subroutine init_field_indices()

    ! tie-gcm
    use fields_module,&
      only: levd0, levd1
    use mpi_module,&
      only: ntask, tasks, mytid
    use params_module,&
      only: glon, glat

    implicit none
    integer :: j

    ! TN and some other fields lowest level is constrainted by lower boundary conditions
    levX0 = levd0
    !ATTENTION level nlevp1 (the highest) can cause problems: nan (e.g., TN) / wrong values(e.g, DEN)
    ! do not use
    levX1 = levd1-1
    nlevX = levX1-levX0+1

    allocate(idx_intern(0:ntask-1))
    allocate(idx_nc(0:ntask-1))
    do j = 0, ntask-1
      call get_lon_idx_intern( j, idx_intern(j)%lon0, idx_intern(j)%lon1, idx_intern(j)%nlons )
      idx_intern(j)%lat0 = tasks(j)%lat0
      idx_intern(j)%lat1 = tasks(j)%lat1
      idx_intern(j)%nlats = tasks(j)%nlats

      idx_nc(j) = idx_intern(j)
      ! index in netcdf file starts at 1
      idx_nc(j)%lon0 =  idx_intern(j)%lon0 - 2;
      idx_nc(j)%lon1 =  idx_intern(j)%lon1 - 2;
    end do

    if(mytid==0) then
      write(*,*) 'intern ---------'
      write(*,*) 'rank lon0 lon1 nlon |  lat0 lat1 nlat'
      do j = 0, ntask-1
        write(*,'(i5,3i5,a,3i5)')  j, idx_intern(j)%lon0, idx_intern(j)%lon1, idx_intern(j)%nlons , &
                              ' | ', idx_intern(j)%lat0, idx_intern(j)%lat1, idx_intern(j)%nlats
      end do
      write(*,*) 'netcdf ---------'
      write(*,*) 'rank lon0 lon1 nlon |  lat0 lat1 nlat |  longitude (deg)   latitude (deg)'
      do j = 0, ntask-1
        write(*,'(i5,3i5,a,3i5,a,4f8.1)')  j, idx_nc(j)%lon0, idx_nc(j)%lon1, idx_nc(j)%nlons , &
                              ' | ', idx_nc(j)%lat0, idx_nc(j)%lat1, idx_nc(j)%nlats , &
                              ' | ', glon(idx_nc(j)%lon0), glon(idx_nc(j)%lon1), &
                                      glat(idx_nc(j)%lat0), glat(idx_nc(j)%lat1)
      end do
    end if

    latX0=idx_intern(mytid)%lat0
    latX1=idx_intern(mytid)%lat1
    lonX0=idx_intern(mytid)%lon0
    lonX1=idx_intern(mytid)%lon1
    nlonX=idx_intern(mytid)%nlons
    nlatX=idx_intern(mytid)%nlats

  end subroutine init_field_indices

  !> Builds this process's local longitude/latitude coordinate vectors (lon_p, lat_p)
  !! from TIE-GCM's global grid.
  subroutine init_horizontal_coords_vector

    ! tie-gcm
    use mpi_module, only: mytid
    use params_module,& 
      only: glon, & ! all longitude coordinates (real) / without periodic cells
            glat     ! all latitude coordinates (real)

    implicit none

    allocate( lon_p, source=glon( idx_nc(mytid)%lon0 : idx_nc(mytid)%lon1 ) )
    allocate( lat_p, source=glat( idx_nc(mytid)%lat0 : idx_nc(mytid)%lat1 ) )

    write(*,'(a,*(f7.1))') "longitudes processed by this rank: ", lon_p
    write(*,'(a,*(f6.1))') "latitudes processed by this rank: ", lat_p

  end subroutine init_horizontal_coords_vector

  !> Deallocates the all allocated variables of the state module.
  subroutine deallocate_state_module

    implicit none

    if (allocated (idx_intern))   deallocate (idx_intern)
    if (allocated (idx_nc))   deallocate (idx_nc)

    if (allocated (lon_p))   deallocate (lon_p)
    if (allocated (lat_p))   deallocate (lat_p)

    call state_vector%destroy
  end subroutine deallocate_state_module

  !> Initializes the mapping between TIE-GCM field names (plus optional calibration
  !! parameters) and their positions in the PDAF augmented state vector.
  !!
  !! Called collectively by each process.
  subroutine state_vector_mapping_init(this, fields, name, calibrations, print_map )

    ! tie-gcm
    use fields_module, &
      only: f4d
    use mpi_module,&
      only: ntask, tasks, mytid

    ! intern
    use model_parameter_handling_module,&
      only: model_parameters
    use uset_module,&
      only: char_uset, str_len_char_uset

    implicit none

    ! args
    class(state_vector_mapping):: this
    character(len=*), dimension(:), intent(in):: fields
    character(len=*), intent(in):: name
    character(len=*), dimension(:), intent(in), optional :: calibrations
    logical, optional, intent(in) :: print_map

    ! local
    integer :: i, j
    integer :: idx

    integer, dimension(:,:), allocatable :: F_R_size
    type(char_uset) :: aug_vec_state, aug_vec_calibration, aug_vec_all
    character(len=str_len_char_uset), allocatable, dimension(:) :: items

    if( size(fields) < 1) then
      call shutdown("construct_state_vector: empty field vector")
    end if

    do i=1, size(fields)
      select case(trim(fields(i)))
        case ("NE","O2P","N2D","TE","TI","OMEGA","POTEN")
          call aug_vec_state%add(item=fields(i))
        case default
          call aug_vec_state%add(item=fields(i))
          call aug_vec_state%add(item=trim(fields(i))//"_NM")
      end select
    end do

    call aug_vec_all%add(other=aug_vec_state)

    if(present(calibrations)) then
      call aug_vec_calibration%add(items=calibrations)
      call aug_vec_all%add(other=aug_vec_calibration)
    end if

    call aug_vec_all%to_array(items)

    allocate( this%fd_idx(1:aug_vec_state%size()) )

    ! find dimensions and idx of fields
    do i=1,aug_vec_state%size()
      ! search in 4d fields
      this%fd_idx(i) = findloc( f4d%short_name, trim(aug_vec_state%at(i)),dim=1 )

      if(  this%fd_idx(i) < 1 ) then
        write(*,'(a,a)') 'Error in init_state. Could not find field: ', aug_vec_state%at(i)
      end if
    end do

    this%idx_f3d_0=1
    this%idx_f3d_1=aug_vec_state%size()

    if(present(calibrations))then
      this%idx_cal_0=aug_vec_state%size()+1
      this%idx_cal_1=aug_vec_all%size()
    else
      this%idx_cal_0=0
      this%idx_cal_1=0
    end if

    ! compute size for tiegcm fields
    allocate( F_R_size(1:aug_vec_all%size(), 0:ntask-1) )
    F_R_size = 0
    do j = 0, ntask-1
      F_R_size(this%idx_f3d_0:this%idx_f3d_1,j) = idx_intern(j)%nlons * tasks(j)%nlats * nlevX
    end do

    ! look up size of calibration parameters
    allocate( this%cal_idx(this%idx_cal_0:this%idx_cal_1) )
    do i=1, aug_vec_calibration%size()
      idx=findloc( model_parameters%name, trim(aug_vec_calibration%at(i)),dim=1 )
      if(  idx < 1 ) then
        write(*,'(a,a)') 'Error in init_state. Could not find calibration parameter: ', aug_vec_calibration%at(i)
      end if
      this%cal_idx(this%idx_cal_0+(i-1))=idx
      F_R_size(this%idx_cal_0+(i-1),0) = model_parameters(idx)%size
    end do

    this%map = construct_mapping(F_R_size, items, name)

    if(present(print_map))then
      if(print_map) then
        if(mytid==0) call this%map%print
      end if
    end if

    call aug_vec_state%to_array(this%tgcm_field_names)
    call aug_vec_calibration%to_array(this%calibration_names)

    ! deallocate local variables
    call aug_vec_state%deallocate()
    call aug_vec_calibration%deallocate()
    call aug_vec_all%deallocate()

    if(allocated(items)) deallocate(items)

    if(allocated(F_R_size)) deallocate(F_R_size)

  end subroutine state_vector_mapping_init

  !> Deallocates a state_vector_mapping instance
  subroutine state_vector_mapping_destroy(this)
    implicit none

   ! arguments
    class(state_vector_mapping), intent(inout) :: this

    if (allocated (this%fd_idx)) deallocate (this%fd_idx)
    call this%map%deallocate()
    if (allocated (this%tgcm_field_names)) deallocate (this%tgcm_field_names)
    if (allocated (this%calibration_names)) deallocate (this%calibration_names)
  end subroutine state_vector_mapping_destroy


  !> Copies this process's subdomain of TIE-GCM field data (at time index itx) into
  !! the PE-local PDAF state vector.
  subroutine state_vector_mapping_fill_state_p(this, state_p, itx)

      ! tie-gcm
      use fields_module,&
        only: f4d
      use mpi_module,&
        only: mytid

      implicit none

      ! args
      class(state_vector_mapping), intent(in) :: this
      real, dimension( : ), contiguous, intent(inout) :: state_p  ! PDAF local state vector
      integer, intent(in) :: itx ! current or previous time step (itc or itp)

      ! local
      integer :: i

      do i = this%idx_f3d_0, this%idx_f3d_1
        state_p( this%map%idx_R(i, mytid)%begin_p : this%map%idx_R(i, mytid)%back_p) =  &
        reshape(  f4d(this%fd_idx(i))%data(        levX0 : levX1,                &
                            idx_intern(mytid)%lon0 : idx_intern(mytid)%lon1, &
                            idx_intern(mytid)%lat0 : idx_intern(mytid)%lat1, &
                                                              itx        ), &
                                        (/ this%map%size(fid=i,rank=mytid)  /)) 
      end do

  end subroutine state_vector_mapping_fill_state_p

  !> Copies TIE-GCM's model-parameter/calibration input constants into
  !! the calibration portion of the state vector, on model rank 0.
  subroutine state_vector_mapping_fill_state_p_calibration(this,state_p)

    ! tie-gcm
    use aurora_module,&
      only: alfac, alfad
    use chemrates_module,&
      only: beta4, beta6, beta7, beta1_0, beta3_0, beta5_0, beta8_0, beta9_0
    use cons_module,&
      only: difk, exp_diff_fac_o2, exp_diff_fac_o1, exp_diff_fac_n2
    use input_module,&
      only: colfac, joulefac
    use mpi_module,&
      only: mytid
    use qrj_module,&
      only: afac

    ! intern
    use array_print_module,&
      only: printMat


    implicit none

    ! arguments
    class(state_vector_mapping), intent(in) :: this
    real, intent(inout), contiguous :: state_p(:)  ! local state vector

    ! local
    integer :: i

    if(mytid==0)then
      do i= this%idx_cal_0, this%idx_cal_1
        SELECT CASE (this%map%fd_name(i))
        CASE ("eddy diffusion")
          state_p(this%map%idx_R(i, 0)%begin_p : &
                  this%map%idx_R(i, 0)%back_p) = difk(:,1)
!           call printMat(state_p(this%map%idx_R(i, 0)%begin_p : &
!                         this%map%idx_R(i, 0)%back_p), "collected eddy" )
!           call printMat(difk(:,1), "difk(:,1)" )
        CASE ("exp_diff_fac_o2")
          state_p(this%map%idx_R(i, 0)%begin_p) = exp_diff_fac_o2
        CASE ("exp_diff_fac_o1")
          state_p(this%map%idx_R(i, 0)%begin_p) = exp_diff_fac_o1
        CASE ("exp_diff_fac_n2")
          state_p(this%map%idx_R(i, 0)%begin_p) = exp_diff_fac_n2
        CASE ("beta1")
          state_p(this%map%idx_R(i, 0)%begin_p) = beta1_0
        CASE ("beta3")
          state_p(this%map%idx_R(i, 0)%begin_p) = beta3_0
        CASE ("beta4")
          state_p(this%map%idx_R(i, 0)%begin_p) = beta4
        CASE ("beta5")
          state_p(this%map%idx_R(i, 0)%begin_p) = beta5_0
        CASE ("beta6")
          state_p(this%map%idx_R(i, 0)%begin_p) = beta6
        CASE ("beta7")
          state_p(this%map%idx_R(i, 0)%begin_p) = beta7
        CASE ("beta8")
          state_p(this%map%idx_R(i, 0)%begin_p) = beta8_0
        CASE ("beta9")
          state_p(this%map%idx_R(i, 0)%begin_p) = beta9_0
        CASE ("colfac")
          state_p(this%map%idx_R(i, 0)%begin_p) = colfac
        CASE ("alfac")
          state_p(this%map%idx_R(i, 0)%begin_p) = alfac
        CASE ("alfad")
          state_p(this%map%idx_R(i, 0)%begin_p) = alfad
        CASE ("joulefac")
          state_p(this%map%idx_R(i, 0)%begin_p) = joulefac
        CASE ("afac","euvafac")
          state_p(this%map%idx_R(i, 0)%begin_p :&
                  this%map%idx_R(i, 0)%back_p) = afac
        CASE DEFAULT
          call shutdown('cannot collect '// this%map%fd_name(i))
        END SELECT
      end do
    end if

  end subroutine state_vector_mapping_fill_state_p_calibration

  !> Writes the (possibly updated/constrained) calibration portion of the state
  !! vector back into TIE-GCM's model-parameter input constants, on model rank 0.
  subroutine state_vector_mapping_distribute_calibration(this,state_p)

    ! tie-gcm
    use aurora_module,&
      only: set_alfac, set_alfad
    use chemrates_module,&
      only: beta2, beta4, beta6, beta7, beta1_0, beta3_0, beta5_0, beta8_0, beta9_0
    use cons_module,&
      only: difk,ndays,exp_diff_fac_o2,exp_diff_fac_o1,exp_diff_fac_n2
    use mpi_module,&
      only: mytid
    use input_module,&
      only: colfac, joulefac
    use qrj_module,&
      only: afac

    ! intern
    use array_print_module,&
      only: printMat
    use model_parameter_handling_module,&
      only: model_parameters

    implicit none

    ! arguments
    class(state_vector_mapping), intent(in) :: this
    real, intent(inout), target, contiguous :: state_p(:)

    ! local
    real, dimension(:), pointer, contiguous :: aug_vector_sub

    integer :: i

    if(mytid==0)then
      do i= this%idx_cal_0, this%idx_cal_1

        aug_vector_sub => state_p(this%map%idx_R(i, 0)%begin_p : &
                                  this%map%idx_R(i, 0)%back_p)

        call model_parameters( this%cal_idx(i) )%constrain(aug_vector_sub)

        SELECT CASE (this%map%fd_name(i))
        CASE ("eddy diffusion")

          difk(:,:) = spread(aug_vector_sub, 2, ndays)
  !         call printMat(state_p(this%map%idx_R(i, 0)%begin_p : &
  !               this%map%idx_R(i, 0)%back_p), "distributed eddy" )
        CASE ("exp_diff_fac_o2")
          exp_diff_fac_o2 = aug_vector_sub(1)
        CASE ("exp_diff_fac_o1")
          exp_diff_fac_o1 = aug_vector_sub(1)
        CASE ("exp_diff_fac_n2")
          exp_diff_fac_n2 = aug_vector_sub(1)
        CASE ("beta1")
          beta1_0 = aug_vector_sub(1)
        CASE ("beta2")
          beta2 = aug_vector_sub(1)
        CASE ("beta3")
          beta3_0 = aug_vector_sub(1)
        CASE ("beta4")
          beta4 = aug_vector_sub(1)
        CASE ("beta5")
          beta5_0 = aug_vector_sub(1)
        CASE ("beta6")
          beta6 = aug_vector_sub(1)
        CASE ("beta7")
          beta7 = aug_vector_sub(1)
        CASE ("beta8")
          beta8_0 = aug_vector_sub(1)
        CASE ("beta9")
          beta9_0 = aug_vector_sub(1)
        CASE ("colfac")
          colfac = aug_vector_sub(1)
        CASE ("alfac")
          call set_alfac(aug_vector_sub(1))
        CASE ("alfad")
          call set_alfad(aug_vector_sub(1))
        CASE ("joulefac")
          joulefac = aug_vector_sub(1)
        CASE ("afac","euvafac")
          afac = aug_vector_sub
        CASE DEFAULT
          call shutdown('cannot distribute '// this%map%fd_name(i))
        END SELECT

      end do

    end if

  end subroutine

  !> Returns a pointer to a named field within the local state vector, reshaped as
  !! a 3D (lev, lon, lat) array.
  function state_vector_mapping_get_field_3d(this,state_p, name) result(field)

    ! tie-gcm
    use mpi_module, only: mytid

    ! intern
    use array_mapping_module, only: map_to_3d

    implicit none

    ! args
    class(state_vector_mapping), intent(in) :: this
    real, dimension(:), target, contiguous, intent(in) :: state_p  ! PDAF local state vector
    character(len=*), intent(in) :: name

    real, dimension(:,:,:), contiguous, pointer:: field

    ! local
    integer :: fid

    fid = findloc( this%map%fd_name, trim(name), dim=1 )
    if(  fid < 1 ) then
      write(*,'(a,a)') name, "is not part of state vector"
      call shutdown(trim(name)//' is not part of state vector')
    end if

    ! ATTENTION assume R order
    call map_to_3d(state_p,field,&
      mat_lb=(/levX0,idx_intern(mytid)%lon0,idx_intern(mytid)%lat0/),&
      mat_ub=(/levX1,idx_intern(mytid)%lon1,idx_intern(mytid)%lat1/),&
      first=this%map%idx_R(fid, mytid)%begin_p,&
      last=this%map%idx_R(fid, mytid)%back_p)

  end function state_vector_mapping_get_field_3d

  !> Returns a flat (1D, unreshaped) pointer to a named field within the local
  !! state vector.
  function state_vector_mapping_get_field(this,state_p, name) result(field)

    ! tie-gcm
    use mpi_module, only: mytid

    ! intern
    use array_mapping_module, only: map_to_3d

    implicit none

    ! args
    class(state_vector_mapping), intent(in) :: this
    real, dimension(:), target, contiguous, intent(in) :: state_p  ! PDAF local state vector
    character(len=*), intent(in) :: name

    real, dimension(:), pointer:: field

    ! local
    integer :: fid

    fid = findloc( this%map%fd_name, trim(name), dim=1 )
    if(  fid < 1 ) then
      write(*,'(a,a)') name, "is not part of state vector"
      call shutdown(trim(name)//' is not part of state vector')
    end if

    ! ATTENTION assume R order
    field => state_p(this%map%idx_R(fid, mytid)%begin_p :&
                      this%map%idx_R(fid, mytid)%back_p)

  end function state_vector_mapping_get_field

!   subroutine get_local_fields(state_p, tn, he, o1, o2, ne,&
!                               tn_nm, he_nm, o1_nm, o2_nm)
! 
!         implicit none
! 
!         ! args
!         real, dimension(:), target, intent(in) :: state_p  ! PDAF local state vector
!         real, dimension(:,:,:), pointer, intent(out), optional :: tn,he,o1,o2,ne,&
!                                                                  tn_nm,he_nm,o1_nm,o2_nm
! 
!         if(present(tn))then
!           tn=>exs(state_p, "TN")
!         end if
!         if(present(he))then
!           he=>exs(state_p, "HE")
!         end if
!         if(present(o1))then
!           o1=>exs(state_p, "O1")
!         end if
!         if(present(o2))then
!           o2=>exs(state_p, "O2")
!         end if
!         if(present(ne))then
!           ne=>exs(state_p, "NE")
!         end if
!         if(present(tn_nm))then
!           tn_nm=>exs(state_p, "TN_NM")
!         end if
!         if(present(he_nm))then
!           he_nm=>exs(state_p, "HE_NM")
!         end if
!         if(present(o1_nm))then
!           o1_nm=>exs(state_p, "O1_NM")
!         end if
!         if(present(o2_nm))then
!           o2_nm=>exs(state_p, "O2_NM")
!         end if
! 
!   end subroutine get_local_fields


  !> Computes the longitude index range (excluding periodic and ghost cells) that a
  !! given task/rank owns in TIE-GCM's f4d fields.
  subroutine get_lon_idx_intern( rank, lon0, lon1, nloni )
    ! rank corresponds to mytid

    ! tie-gcm
    use mpi_module, only: tasks, ntaski

    implicit none

    ! arguments
    integer, intent(in)  :: rank
    integer, intent(out) :: lon0, lon1, nloni

    ! consider periodic design of grid in longitude
    if (ntaski == 1 ) then
        !  grid is not decomposed in longitude 
        lon0 = tasks(rank)%lon0 +2
        lon1 = tasks(rank)%lon1 -2
    else if( tasks(rank)%mytidi == 0) then
        ! I am at west end (mytidi==0) of task row 
        lon0 = tasks(rank)%lon0 + 2
        lon1 = tasks(rank)%lon1
    else if(  tasks(rank)%mytidi == ntaski-1   ) then
        !  I am at  east end (mytidi==ntaski-1)
        lon0 = tasks(rank)%lon0
        lon1 = tasks(rank)%lon1 -2
    else
        lon0 = tasks(rank)%lon0
        lon1 = tasks(rank)%lon1
    end if

    nloni = lon1-lon0+1

  end subroutine get_lon_idx_intern

  !> Computes the mass-fraction-to-number-density conversion factor for a given
  !! species (O1/O2/N2/HE) from xnmbar.
  subroutine calc_conversion_factors_mmr_to_n(xnmbar, species, factor)

    ! tie-gcm
    use cons_module,&
      only: rmassinv_o1, rmassinv_o2, rmassinv_n2, rmassinv_he
    use fields_module,&
      only: levd0,levd1,lond0,lond1,latd0,latd1

    implicit none

    ! arguments
    real, dimension(levd0:levd1,lond0:lond1,latd0:latd1), intent(in):: xnmbar
    character(len=*), intent(in) :: species
    real, dimension(levd0:levd1,lond0:lond1,latd0:latd1), intent(out):: factor

    ! local
    real :: molarmass_inv
    logical :: apply

    apply = .true.
    select case ( species )
        case ('O1', 'O1_NM')
            molarmass_inv = rmassinv_o1
        case ('O2', 'O2_NM')
            molarmass_inv = rmassinv_o2
        case ('N2', 'N2_NM')
            molarmass_inv = rmassinv_n2
        case ('HE', 'HE_NM')
            molarmass_inv = rmassinv_he
        case default
            write(*,*) 'error species ', species, ' can not be converted from mass fraction to number density'
            apply = .false.
    end select

    if(apply) then
      factor = xnmbar * molarmass_inv
    else
      factor = 1.0
    end if

  end subroutine calc_conversion_factors_mmr_to_n

end module state_module
