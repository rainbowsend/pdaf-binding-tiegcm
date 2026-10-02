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
! Registry of all computable/output quantities (TIE-GCM fields, densities, gradients, geometric height) with routines to compute each from TIE-GCM or state-vector data.

module quantity_info_module

  ! tie-gcm
  use fields_module, only: longname_len, shortname_len, units_len

  implicit none

  integer, parameter :: LEVEL_INT = 0  ! on interfaces
  integer, parameter :: LEVEL_MID = 1  ! on midpoints
  integer, parameter :: LEVEL_NONE = 2 ! no vertical extend, e.g., integrated quantities

  integer, parameter :: STEP_CURRENT = 0
  integer, parameter :: STEP_PREVIOUS = 1

  ! 1 tecu = 1e+16 electrons/m2 = 1e+12 electrons/cm2
  real, parameter :: cm2_per_tecu = 1e+12


  ! ATTENTION data of quantity_info is either stored in a f4d TIE-GCM field or
  ! in member variable additional_fields (initalized in init_pdaf)
  type :: quantity_info
    character(len=shortname_len) :: name
    character(len=units_len) :: unit
    character(len=longname_len) :: long_name
    integer :: level
    integer :: step
    real, dimension(:,:,:), contiguous, pointer :: data => null()
    procedure(calc_interface_f4d), pointer, private :: calc_from_f4d
    procedure(calc_interface_state), pointer, private :: calc_from_state
    contains
    procedure, pass(this), public :: level_string => quantity_info_level_string
    procedure, pass(this), public :: calc => quantity_info_calc
  end type

  type :: quantity_type
    type(quantity_info), pointer :: info
  end type

  type(quantity_info), dimension(:), allocatable, target :: quantity_infos

  real, dimension(:,:,:,:), allocatable, target :: additional_fields ! allocated in init_pdaf
  ! Quantities without vertical extent (LEVEL_NONE) are stored separately with
  ! a single level.
  ! Their data pointer keeps rank three, thus they are used just like any
  ! other quantity.
  real, dimension(:,:,:,:), allocatable, target :: additional_fields_2d ! allocated in init_pdaf

  abstract interface
    !> Abstract interface for computing a quantity's data from the current TIE-GCM fields.
    subroutine calc_interface_f4d(this)
      import quantity_info
      class(quantity_info) :: this
    end subroutine
    !> Abstract interface for computing a quantity's data from the PDAF state vector.
    subroutine calc_interface_state(this, state_map, state_p)
      use state_module, only : state_vector_mapping
      import quantity_info
      class(quantity_info) :: this
      type(state_vector_mapping), intent(in) :: state_map
      real, dimension(:), contiguous, intent(in) :: state_p
    end subroutine
  end interface

  contains

  !> Points this%data at the current-timestep slice of the matching TIE-GCM f4d field.
  subroutine get_from_f4d(this)

    ! tiegcm
    use fields_module, only: f4d,itc

    implicit none

    class(quantity_info) :: this

    integer :: fidx

    ! state has same quantity as observation
    fidx = findloc(f4d%short_name,this%name,dim=1)

    if(fidx<lbound(f4d,dim=1)) then
      call shutdown(this%name//' is not a tiegcm field')
    end if

    this%data => f4d(fidx)%data(:,:,:,itc)

  end subroutine

  !> Fills this%data from the matching field in the PDAF state vector, reshaped to the local 3D grid.
  subroutine get_from_state(this, state_map, state_p)

    ! intern
    use state_module, only : state_vector_mapping, levX0, levX1, lonX0, lonX1, latX0, latX1, nlevX, nlonX, nlatX

    implicit none

    class(quantity_info) :: this
    type(state_vector_mapping), intent(in) :: state_map
    real, dimension(:), contiguous, intent(in) :: state_p

    this%data( levX0 : levX1, &
               lonX0 : lonX1, &
               latX0 : latX1) &
        = reshape( state_map%field(state_p, this%name), &
                  (/nlevX, nlonX, nlatX /))

  end subroutine


  !> Computes neutral mass density at midpoints from the current-timestep TIE-GCM tn/o2/o1/he fields.
  subroutine calc_den(this)

    ! tiegcm
    use fields_module, only: levd0,levd1,lond0,lond1,latd0,latd1,tn,o2,o1,he,itc
    use aerostatic_diag, only: calc_diag

    implicit none

    class(quantity_info) :: this

    call calc_diag(lon0=lond0,&
                   lon1=lond1,&
                   lev0=levd0,&
                   lev1=levd1,&
                   lat0=latd0,&
                   lat1=latd1,&
                   tn=tn(:,:,:,itc),&
                   o2=o2(:,:,:,itc),&
                   o1=o1(:,:,:,itc),&
                   he=he(:,:,:,itc),&
                   rhomid=this%data)

  end subroutine

  !> Computes neutral mass density at midpoints from the previous-timestep (_NM) TIE-GCM fields.
  subroutine calc_den_nm(this)

    ! tiegcm
    use fields_module, only: levd0,levd1,lond0,lond1,latd0,latd1,tn_nm,o2_nm,o1_nm,he_nm,itc
    use aerostatic_diag, only: calc_diag

    implicit none

    class(quantity_info) :: this

    call calc_diag(lon0=lond0,&
                   lon1=lond1,&
                   lev0=levd0,&
                   lev1=levd1,&
                   lat0=latd0,&
                   lat1=latd1,&
                   tn=tn_nm(:,:,:,itc),&
                   o2=o2_nm(:,:,:,itc),&
                   o1=o1_nm(:,:,:,itc),&
                   he=he_nm(:,:,:,itc),&
                   rhomid=this%data)

  end subroutine

  !> Computes neutral mass density at midpoints from the current-timestep fields stored in the PDAF state vector.
  subroutine calc_den_from_state(this, state_map, state_p)

    use state_module, only : state_vector_mapping, levX0, levX1, lonX0, lonX1, latX0, latX1
    use aerostatic_diag, only: calc_diag

    implicit none

    class(quantity_info) :: this
    type(state_vector_mapping), intent(in) :: state_map
    real, dimension(:), contiguous, intent(in) :: state_p

    call calc_diag(lon0=lonX0,&
                   lon1=lonX1,&
                   lev0=levX0,&
                   lev1=levX1,&
                   lat0=latX0,&
                   lat1=latX1,&
                   tn=state_map%field_3d(state_p,"TN"),&
                   o2=state_map%field_3d(state_p,"O2"),&
                   o1=state_map%field_3d(state_p,"O1"),&
                   he=state_map%field_3d(state_p,"HE"),&
                   rhomid=this%data(levX0:levX1,lonX0:lonX1,latX0:latX1))
  end subroutine

  !> Computes neutral mass density at midpoints from the previous-timestep (_NM) fields stored in the PDAF state vector.
  subroutine calc_den_nm_from_state(this, state_map, state_p)

    use state_module, only : state_vector_mapping, levX0, levX1, lonX0, lonX1, latX0, latX1
    use aerostatic_diag, only: calc_diag

    implicit none

    class(quantity_info) :: this
    type(state_vector_mapping), intent(in) :: state_map
    real, dimension(:), contiguous, intent(in) :: state_p

    call calc_diag(lon0=lonX0,&
                   lon1=lonX1,&
                   lev0=levX0,&
                   lev1=levX1,&
                   lat0=latX0,&
                   lat1=latX1,&
                   tn=state_map%field_3d(state_p,"TN_NM"),&
                   o2=state_map%field_3d(state_p,"O2_NM"),&
                   o1=state_map%field_3d(state_p,"O1_NM"),&
                   he=state_map%field_3d(state_p,"HE_NM"),&
                   rhomid=this%data(levX0:levX1,lonX0:lonX1,latX0:latX1))
  end subroutine

!   subroutine calc_zg_mid_from_state(this, state_map, state_p)
!     use state_module, only : state_vector_mapping, levX0, levX1, lonX0, lonX1, latX0, latX1
!     use aerostatic_diag, only: calc_diag, int_to_mid
!
!     implicit none
!
!     class(quantity_info) :: this
!     type(state_vector_mapping), intent(in) :: state_map
!     real, dimension(:), contiguous, intent(in) :: state_p
!
!     call calc_diag(lon0=lonX0,&
!                    lon1=lonX1,&
!                    lev0=levX0,&
!                    lev1=levX1,&
!                    lat0=latX0,&
!                    lat1=latX1,&
!                    tn=state_map%field_3d(state_p,"TN"),&
!                    o2=state_map%field_3d(state_p,"O2"),&
!                    o1=state_map%field_3d(state_p,"O1"),&
!                    he=state_map%field_3d(state_p,"HE"),&
!                    zg=zg_mid(levX0:levX1,lonX0:lonX1,latX0:latX1))
!     zg_mid(levX0:levX1,lonX0:lonX1,latX0:latX1) = &
!       int_to_mid(lonX0,lonX1,levX0,levX1,latX0,latX1,&
!       zg_mid(levX0:levX1,lonX0:lonX1,latX0:latX1))
!   end subroutine

!   subroutine calc_ne(this)
!     ! tiegcm
!     use fields_module, only: levd0,levd1,lond0,lond1,latd0,latd1,op,o2p,nop,nplus,n2p,itc
!
!     implicit none
!
!     class(quantity_info) :: this
!
!     ! this overrides ne in fields module
!     this%data = op(:,:,:,itc)+o2p(:,:,:,itc)+nop(:,:,:)+nplus(:,:,:)+n2p(:,:,:)
!   end subroutine

  !> Integrates the electron density over the vertical column.
  !!
  !! Both the electron density and the geometric height are located at
  !! interfaces, thus the trapezoidal rule is applied between them.
  !!
  !! ATTENTION the integral only covers the levels of TIE-GCM, which reaches
  !! from about 90 km up to 500-700 km, depending on ZITOP and solar activity.
  !! The plasmasphere above the upper boundary is not represented, although it
  !! contributes a considerable part (typically 10-30 %) of the total electron
  !! content observed by GNSS. Electrons below the lower boundary are missing
  !! as well, but their contribution is negligible in comparison. Consequently
  !! the result is biased low with respect to such observations.
  function electron_content(ne, zg) result(column)

    implicit none

    ! arguments
    real, dimension(:,:,:), intent(in) :: ne ! electron density at interfaces (1/cm3)
    real, dimension(:,:,:), intent(in) :: zg ! geometric height at interfaces (cm)

    ! result
    real, dimension(size(ne,dim=2),size(ne,dim=3)) :: column ! electron content (1/cm2)

    ! local
    integer :: k

    ! TIE-GCM uses a staggered vertical grid. Both fields integrated here are
    ! located at the interfaces, which bound the layers:
    !   ne  electron density (1/cm3), vcoord "interfaces" (see fields.F)
    !   zg  geometric height (cm)   , at interfaces too   (see addiag.F)
    !
    ! This is exactly what the trapezoidal rule needs, since integrating over
    ! a single layer requires its thickness and the value enclosed by it:
    !
    !   interface k+1  ---------  ne(k+1), zg(k+1)   ^
    !   midpoint       - - - - -                     |  dz = zg(k+1)-zg(k)
    !   interface k    ---------  ne(k)  , zg(k)     v
    !
    ! The thickness follows from the difference of the two bounding
    ! interfaces. The mean of the two bounding interface values is the value
    ! at the enclosed midpoint.

    column = 0.
    do k = 1, size(ne,dim=1)-1
      column = column + 0.5*(ne(k,:,:)+ne(k+1,:,:)) * (zg(k+1,:,:)-zg(k,:,:))
    end do

  end function electron_content

  !> Computes the vertical total electron content from the current-timestep
  !! TIE-GCM electron density and geometric height.
  subroutine calc_vtec(this)

    ! tiegcm
    use fields_module, only: levd0,ne,zg,itc

    ! intern
    use state_module, only: levX0, levX1

    implicit none

    class(quantity_info) :: this

    ! the uppermost TIE-GCM level holds invalid values, thus levX1 is the
    ! highest level that can be integrated
    this%data(levd0,:,:) = electron_content(ne(levX0:levX1,:,:,itc),&
                                            zg(levX0:levX1,:,:)) / cm2_per_tecu

  end subroutine

  !> Computes the vertical total electron content from the electron density in
  !! the PDAF state vector and the geometric height derived from it.
  subroutine calc_vtec_from_state(this, state_map, state_p)

    ! tiegcm
    use fields_module, only: levd0

    ! intern
    use aerostatic_diag, only: calc_diag
    use state_module, only : state_vector_mapping, levX0, levX1, lonX0, lonX1, latX0, latX1

    implicit none

    class(quantity_info) :: this
    type(state_vector_mapping), intent(in) :: state_map
    real, dimension(:), contiguous, intent(in) :: state_p

    ! local
    real, dimension(levX0:levX1,lonX0:lonX1,latX0:latX1) :: zg_state

    ! the geometric height depends on the assimilated fields, thus it has to be
    ! derived from the state vector instead of using the one of the model
    call calc_diag(lon0=lonX0,&
                   lon1=lonX1,&
                   lev0=levX0,&
                   lev1=levX1,&
                   lat0=latX0,&
                   lat1=latX1,&
                   tn=state_map%field_3d(state_p,"TN"),&
                   o2=state_map%field_3d(state_p,"O2"),&
                   o1=state_map%field_3d(state_p,"O1"),&
                   he=state_map%field_3d(state_p,"HE"),&
                   zg=zg_state)

    this%data(levd0,lonX0:lonX1,latX0:latX1) = &
      electron_content(state_map%field_3d(state_p,"NE"), zg_state) / cm2_per_tecu

  end subroutine

  !> Computes molecular nitrogen mass mixing ratio from the current-timestep TIE-GCM o2/o1/he fields.
  subroutine calc_n2(this)

    ! tiegcm
    use fields_module, only: o2,o1,he,itc
    use aerostatic_diag, only: molecular_nitrogen

    implicit none

    class(quantity_info) :: this

    this%data = molecular_nitrogen(o2(:,:,:,itc), o1(:,:,:,itc), he(:,:,:,itc))
  end subroutine

  !> Computes a quantity's time gradient from its current and previous-timestep TIE-GCM f4d fields.
  subroutine calc_grad(this)

    ! tiegcm
    use fields_module, only: f4d, itc
    use aerostatic_diag, only: calc_diag

    implicit none

    class(quantity_info) :: this

    integer :: charlen
    integer :: idx, idx_nm

    charlen = len_trim(this%name)

    idx = findloc(f4d%short_name,this%name(1:charlen-2),dim=1)
    idx_nm = findloc(f4d%short_name,this%name(1:charlen-2)//'_NM',dim=1)

    this%data = calc_gradient(f4d(idx)%data(:,:,:,itc), f4d(idx_nm)%data(:,:,:,itc))

  end subroutine

  !> Computes a quantity's time gradient from its current and previous-timestep values stored in the PDAF state vector.
  subroutine calc_grad_from_state(this, state_map, state_p)

    ! intern
    use state_module, only : state_vector_mapping, levX0, levX1, lonX0, lonX1, latX0, latX1

    implicit none

    class(quantity_info) :: this
    type(state_vector_mapping), intent(in) :: state_map
    real, dimension(:), contiguous, intent(in) :: state_p

    ! local
    integer :: charlen

    real, dimension(:,:,:), pointer, contiguous :: current, previous

    charlen = len_trim(this%name)

    current => state_map%field_3d(state_p,this%name(1:charlen-2))
    previous => state_map%field_3d(state_p,this%name(1:charlen-2)//"_NM")

    this%data( levX0 : levX1, &
               lonX0 : lonX1, &
               latX0 : latX1) &
        =  calc_gradient(current, previous)

  end subroutine

  !> True if the quantity can be interpolated to a single point.
  !!
  !! Currently, quantities without vertical extent (LEVEL_NONE), e.g. VTEC,
  !! cannot, since the sparse interpolator builds a vertical spline per
  !! source column. This is the only definition of that restriction.
  function supports_point_interpolation(varname) result(supported)

    implicit none

    ! arguments
    character(len=*), intent(in) :: varname

    ! result
    logical :: supported

    ! local
    type(quantity_info), pointer :: info

    info => get_info(varname)

    supported = info%level /= LEVEL_NONE

  end function supports_point_interpolation

  !> Looks up and returns a pointer to the quantity_info entry with the given name, aborting if unknown.
  function get_info(varname) result(info)

    implicit none

    ! arguments
    character(len=*), intent(in) :: varname

    ! result
    type(quantity_info), pointer :: info

    ! local
    integer :: idx

    idx = findloc( quantity_infos%name, trim(varname), dim=1 )
    if( idx < 1) then
      call shutdown('unknown quantity: "'//varname//'"')
    end if

    info => quantity_infos(idx)

  end function

  !> Scans a list of quantities and flags which geometric-height variants (mid/int, current/previous) are needed to compute them.
  subroutine get_required_zg(quantities,  req_mid, req_mid_nm, req_int, req_int_nm)

    implicit none

    ! arguments
    type(quantity_type), dimension(:), intent(in) :: quantities
    logical, intent(out) :: req_mid, req_mid_nm, req_int, req_int_nm

    ! local
    integer :: i

    req_mid = .false.
    req_mid_nm = .false.
    req_int = .false.
    req_int_nm = .false.

    do i = 1, size(quantities)
      select case(quantities(i)%info%level)
        case(LEVEL_MID)
          select case(quantities(i)%info%step)
            case(STEP_CURRENT)
              req_mid = .true.
            case(STEP_PREVIOUS)
              req_mid_nm = .true.
          end select
        case(LEVEL_INT)
          select case(quantities(i)%info%step)
            case(STEP_CURRENT)
              req_int = .true.
            case(STEP_PREVIOUS)
              req_int_nm = .true.
          end select
        case(LEVEL_NONE)
          ! quantities without vertical extent are interpolated horizontally
          ! only, thus they do not require any geometric height
      end select
    end do

  end subroutine

  !> Builds the registry of all computable quantities (TIE-GCM prognostic fields, their gradients, density, and geometric height) and wires up each one's compute routines.
  subroutine init_quantity_info_module

    ! tie-gcm
    use fields_module, only: f4d, nf4d

    ! intern
    use character_routines_module, only: string_ends_with

    implicit none

    ! local
    integer :: i,j
    integer :: idx
    character(len=128) :: msg

    ! +11 for gradients of f4d fields
    ! +6 for addional fields that are not stored in f4d
    allocate(quantity_infos(nf4d+11+6))

    i = 0
    do j=1,nf4d
      i = i+1
      write(quantity_infos(i)%name,'(a)') f4d(j)%short_name
      write(quantity_infos(i)%unit,'(a)') f4d(j)%units
      write(quantity_infos(i)%long_name,'(a)') f4d(j)%long_name

      select case (f4d(j)%vcoord)
        case ("midpoints")
          quantity_infos(i)%level=LEVEL_MID
        case ("interfaces")
          quantity_infos(i)%level=LEVEL_INT
      end select

      if(string_ends_with(trim(quantity_infos(i)%name),'_NM'))then
        quantity_infos(i)%step=STEP_PREVIOUS
      else
        quantity_infos(i)%step=STEP_CURRENT
      end if

      quantity_infos(i)%calc_from_f4d => get_from_f4d
      quantity_infos(i)%calc_from_state => get_from_state
    end do

    ! gradients of prognostic fields
    do j=1,nf4d
       if(string_ends_with(trim(f4d(j)%short_name),'_NM').eqv..false.)then
          idx = findloc(f4d%short_name,trim(f4d(j)%short_name)//'_NM',dim=1)
          if(idx > 0) then
            i=i+1
            write(quantity_infos(i)%name,'(a)') trim(f4d(j)%short_name)//"_G"
            write(quantity_infos(i)%unit,'(a)') trim(f4d(j)%units)//" 1/s"
            write(quantity_infos(i)%long_name,'(a)') trim(f4d(j)%long_name)//" gradient"
            select case (f4d(j)%vcoord)
              case ("midpoints")
                quantity_infos(i)%level=LEVEL_MID
              case ("interfaces")
                quantity_infos(i)%level=LEVEL_INT
            end select
            quantity_infos(i)%step=STEP_CURRENT

            quantity_infos(i)%calc_from_f4d => calc_grad
            quantity_infos(i)%calc_from_state => calc_grad_from_state

          end if
       end if
    end do

    i=i+1
    write(quantity_infos(i)%name,'(a)') "DEN"
    write(quantity_infos(i)%unit,'(a)') "g/cm3"
    write(quantity_infos(i)%long_name,'(a)') "neutral mass density"
    quantity_infos(i)%step=STEP_CURRENT
    quantity_infos(i)%level=LEVEL_MID
    quantity_infos(i)%calc_from_f4d => calc_den
    quantity_infos(i)%calc_from_state => calc_den_from_state

    i=i+1
    write(quantity_infos(i)%name,'(a)') "DEN_NM"
    write(quantity_infos(i)%unit,'(a)') "g/cm3"
    write(quantity_infos(i)%long_name,'(a)') "neutral mass density"
    quantity_infos(i)%step=STEP_PREVIOUS
    quantity_infos(i)%level=LEVEL_MID
    quantity_infos(i)%calc_from_f4d => calc_den_nm
    quantity_infos(i)%calc_from_state => calc_den_nm_from_state

    i=i+1
    write(quantity_infos(i)%name,'(a)') "ZG"
    write(quantity_infos(i)%unit,'(a)') "cm"
    write(quantity_infos(i)%long_name,'(a)') "geometric height at interfaces"
    quantity_infos(i)%step=STEP_CURRENT
    quantity_infos(i)%level=LEVEL_INT

    i=i+1
    write(quantity_infos(i)%name,'(a)') "ZGMID"
    write(quantity_infos(i)%unit,'(a)') "cm"
    write(quantity_infos(i)%long_name,'(a)') "geometric height at midpoints"
    quantity_infos(i)%step=STEP_CURRENT
    quantity_infos(i)%level=LEVEL_MID
    quantity_infos(i)%calc_from_f4d => null()
    quantity_infos(i)%calc_from_state => null()

!     i=i+1
!     write(quantity_infos(i)%name,'(a)') "P"
!     write(quantity_infos(i)%unit,'(a)') "bar"
!     write(quantity_infos(i)%long_name,'(a)') "pressure"
!     quantity_infos(i)%step=STEP_CURRENT
!     quantity_infos(i)%level=LEVEL_MID
!     ! TODO implement calc_from_

    i=i+1
    write(quantity_infos(i)%name,'(a)') "N2"
    write(quantity_infos(i)%unit,'(a)') "mmr"
    write(quantity_infos(i)%long_name,'(a)') "molecular nitrogen"
    quantity_infos(i)%step=STEP_CURRENT
    quantity_infos(i)%level=LEVEL_MID
    quantity_infos(i)%calc_from_f4d => calc_n2
    quantity_infos(i)%calc_from_state => null()

    i=i+1
    write(quantity_infos(i)%name,'(a)') "VTEC"
    write(quantity_infos(i)%unit,'(a)') "tecu"
    write(quantity_infos(i)%long_name,'(a)') "vertical total electron content"
    quantity_infos(i)%step=STEP_CURRENT
    quantity_infos(i)%level=LEVEL_NONE
    quantity_infos(i)%calc_from_f4d => calc_vtec
    quantity_infos(i)%calc_from_state => calc_vtec_from_state

    if(i/=size(quantity_infos,dim=1))then
      write(msg,*) 'initalized ', i , ' quantities but expected ', size(quantity_infos,dim=1)
      call shutdown("quantity_infos was not allocated correctly! "//msg)
    end if

    call write_quantities

  end subroutine

  !> Deallocates the quantity registry and the additional-fields array.
  subroutine finalize_quantity_info_module
    implicit none
    if(allocated(quantity_infos)) deallocate(quantity_infos)
    if(allocated(additional_fields)) deallocate(additional_fields)
    if(allocated(additional_fields_2d)) deallocate(additional_fields_2d)
  end subroutine

  !> Computes this quantity's data, from the state vector if given, otherwise from the current TIE-GCM fields.
  subroutine quantity_info_calc(this,state_map,state_p)

    use state_module, only : state_vector_mapping

    implicit none

    ! input arguments
    class(quantity_info), intent(in) :: this
    type(state_vector_mapping), intent(in), optional :: state_map
    real, dimension(:), contiguous, intent(in), optional :: state_p

    if((present(state_map)).and.(present(state_p))) then
      if(associated(this%calc_from_state))then
        call this%calc_from_state(state_map,state_p)
      else
        write(*,*) "WARNING: no operator to calculate ", this%name, " from state vector"
      end if
    else
      if(associated(this%calc_from_f4d))then
        call this%calc_from_f4d()
      else
        write(*,*) "WARNING: no operator to calculate ", this%name, " from tiegcm fields"
      end if
    end if
  end subroutine

  !> Returns "midpoints" or "interfaces" for this quantity's vertical level.
  function quantity_info_level_string(this) result(level)

    implicit none

    ! input arguments
    class(quantity_info), intent(in) :: this

    ! result
    character(len=16) :: level

    select case(this%level)
      case(LEVEL_MID)
      level="midpoints"
      case(LEVEL_INT)
        level="interfaces"
      case(LEVEL_NONE)
        level="none"
      case default
        call shutdown('quantity_info_level_string: unhandled level of quantity '//trim(this%name))
    end select

  end function quantity_info_level_string

  !> Computes the time gradient between two field snapshots by finite-differencing over one TIE-GCM time step.
  function calc_gradient(current, previous) result(grad)

    use input_module, only: step ! length per step in seconds

    implicit none

    ! arguments
    real, dimension(:,:,:), intent(in) :: current
    real, dimension(:,:,:), intent(in) :: previous

    ! result
    real, dimension(size(current,dim=1), size(current,dim=2), size(current,dim=3)) :: grad

    grad = (current - previous)/step

  end function

  !> Prints a table of all registered quantities with their name, long name, and unit.
  subroutine write_quantities
    implicit none

    integer :: i
    write(*,*) 'Table of available quantities that can be written by result file writer:'
    do i = 1, size(quantity_infos,dim=1)
      write(*,'(a16,a80,a16)') quantity_infos(i)%name, quantity_infos(i)%long_name, quantity_infos(i)%unit
    end do
  end subroutine

end module
