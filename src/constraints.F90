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
! Enforces physical bounds/consistency on the analysis state (mass fractions, positivity, temperature order, charge neutrality, NaNs).

module constraints

! intern
use configuration, only: cfg_log

implicit none

contains

!> Applies the configured set of physical constraints (species, temperature, electron density, neutrality, gradients, NaNs) to the current analysis state.
subroutine constrain_state

  use configuration, only: cfg_state, cfg_constraints
  use state_module, only: state_vector

  implicit none

  logical :: major
  logical :: minor

  major = any((/cfg_state%o2,&
                cfg_state%o1,&
                cfg_state%he/))
  minor = any((/cfg_state%atomic_argon,&
                cfg_state%nitric_oxide,&
                cfg_state%excited_atomic_nitrogen_4s,&
                cfg_state%excited_atomic_nitrogen_2d/))

  if( cfg_constraints%relative_gradient_limit > 0 ) then
    call gradient_constraint(state_vector%tgcm_field_names,&
                             cfg_constraints%relative_gradient_limit)
  end if

  if( major .or. minor) then
    call species_constraint(minor)
  end if

  if(cfg_state%tn .or. &
     cfg_state%electron_temperature .or. &
     cfg_state%ion_temperature ) then
    call temperature_constraint()
  end if

  if(cfg_state%ne) then
    call electron_density_constraint
  end if

  if(cfg_state%atomic_oxygen_ion_density) then
    call atomic_oxygen_ion_density_constraint
  end if

  if(cfg_state%molecular_oxygen_ion_density) then
    call molecular_oxygen_ion_density_constraint
  end if

  ! This replaces NE so that it is equal to the sum of all ions
  ! Since only O1+ and O2+ ions are progonsitc one can ignore this constraint, since it will be compensated by the other ions not in the state vector.
  if(cfg_constraints%quasi_neutral_ionosphere) then
    if(cfg_state%molecular_oxygen_ion_density .or. &
      cfg_state%atomic_oxygen_ion_density ) then

      if(cfg_state%ne) then
        write(*,*) "!!! WARNING !!!  NE is in state vector. But it is overridden by neutrallity constriant ", &
                "NE = sum(ions)"
      end if

      call electrically_neutral_constraint(ions_are_ref=.true.)
    end if
  end if

  call no_nan_constraint

end subroutine

!> Applies the mass-fraction constraint to major (and optionally minor) species, for both the previous and current time levels.
subroutine species_constraint(include_minor)

  ! extern
  use mpi_f08

  ! tie-gcm
  use fields_module,&
      only: i_o1,i_o2,i_he,i_o1_nm,i_o2_nm,i_he_nm, n2d, itc
  use fields_module,&
    only: i_ar, i_n4s, i_no, i_n2d, i_ar_nm, i_n4s_nm, i_no_nm
  implicit none


  logical, intent(in) :: include_minor

  real, dimension(:,:,:), allocatable :: tmp

  if(cfg_log%verbose_level>0) write(*,*) 'apply mass_fraction_constraint for previous epoch'
  if(include_minor) then
    allocate(tmp,source=n2d(:,:,:,itc))
    call mass_fraction_constraint((/i_o1_nm,i_o2_nm,i_he_nm, i_ar_nm, i_n4s_nm, i_no_nm, i_n2d/))
  else
    call mass_fraction_constraint((/i_o1_nm,i_o2_nm,i_he_nm/))
  end if

  if(cfg_log%verbose_level>0) write(*,*) 'apply mass_fraction_constraint for current epoch'
  if(include_minor) then
    n2d(:,:,:,itc) = tmp
    call mass_fraction_constraint((/i_o1,i_o2,i_he, i_ar, i_n4s, i_no, i_n2d/))
  else
    call mass_fraction_constraint((/i_o1,i_o2,i_he/))
  end if

  if(allocated(tmp)) deallocate(tmp)


end subroutine

!> Rescales the given species so their mass fractions sum to at most one in every grid cell.
subroutine mass_fraction_constraint(fids)

  ! tie-gcm
  use fields_module,&
      only: itc, f4d, levd0, levd1, lond0, lond1, latd0, latd1
  use mpi_module,&
      only: mytid
  use params_module,&
      only: nlat, nlon, nlev

  ! intern
  use mod_parallel_pdaf, only: task_id

  implicit none

  ! args
  integer, dimension(:), intent(in) :: fids

  ! local
  real, dimension(levd0:levd1, lond0:lond1, latd0:latd1) :: C
  logical, dimension(levd0:levd1, lond0:lond1, latd0:latd1) :: mask
  integer :: i

  integer :: number_cells_sum_mass_fraction_condition_failed

  real, dimension(:,:,:), pointer, contiguous :: species

  ! To  make sure mass fraction of N2 is larger than 0 the sum of
  ! O1, O2 and He must be less than one.
  ! Thus, substract small value from one
  real, parameter ::  Cmax = 1.0-1E-12

  C = 0
  do i=1, size(fids,dim=1)
    if(cfg_log%verbose_level>0) write(*,*) 'constraint ', f4d(fids(i))%long_name
    species => f4d(fids(i))%data(:,:,:,itc)
    ! mass fraction has to be in (0,1]
    ! to avoid division by zero limit values to a small value but not zero
    ! (in TIE-GCM species occur in denominator of some fractions)
    call limit_values(species,1E-12,1.)

    C = C + species

  end do

  ! Sum of all species has to be smaller than 1

  mask = .false.
  where( C > Cmax) mask = .true.
  do i=1, size(fids,dim=1)
    species => f4d(fids(i))%data(:,:,:,itc)
    species = MERGE(species*Cmax/C, species, mask)
  end do

  call count_failed_cells(mask, number_cells_sum_mass_fraction_condition_failed)

  if( mytid .eq. 0 ) then
    if (number_cells_sum_mass_fraction_condition_failed .gt. 0) then
    write(*,*) 'ensemble member ', task_id, ' WARNING species constraint: rescaling mass fractions in ',&
                number_cells_sum_mass_fraction_condition_failed, '/', nlat*nlon*nlev ,' cells'
    end if

  end if

end subroutine

!> Clamps electron density to a physically plausible range.
subroutine electron_density_constraint()

  ! tie-gcm
  use fields_module, only: ne, itc

  implicit none

   if(cfg_log%verbose_level>0) write(*,*) 'constraint electron density'

  ! see elden.F: "Insure positive Ne (at least 3100):"
  ! max: 1.2 * max value in TIE-GCM hallooween storm run rounded to two significant figures
   call limit_values(ne(:,:,:,itc), 3100., 4E+6)

end subroutine

!> Clamps neutral/electron/ion temperatures to plausible ranges and enforces Te >= Tn and Ti >= Tn.
subroutine temperature_constraint()

  ! tie-gcm
  use fields_module, only: itc, tn, tn_nm, te, ti
  use configuration, only: cfg_state

  implicit none

  if(cfg_log%verbose_level>0) write(*,*) 'constraint temperatures'

  if(cfg_state%tn) then
    ! Temperature in Kelvin, can not be negative
    !  dt.F:377:! Tn must be at least 100 deg:
    !
    ! max value in TIE-GCM data prim files 1582
    ! max n nrlmsis run 2400*1.2
    call limit_values(tn(:,:,:,itc), 100., 2800.)
    call limit_values(tn_nm(:,:,:,itc), 100., 2800.)
  end if

  if(cfg_state%electron_temperature) then
    ! min and max value in TIE-GCM data prim files 142.076 5098.88
    ! max: 1.2 * max value in TIE-GCM hallooween storm run rounded to two significant figures
    call limit_values(te(:,:,:,itc), 100., 6100.)
  end if

  if(cfg_state%ion_temperature) then
    ! min and max value in TIE-GCM data prim files 142.033 3695.77
    ! max: 1.2 * max value in TIE-GCM hallooween storm run rounded to two significant figures
    call limit_values(ti(:,:,:,itc), 100., 4200.)
  end if

  ! settei.F
  ! Te must be >= Tn
  ! ti must be at least as large as tn:

  where(te(:,:,:,itc) < tn(:,:,:,itc)) te(:,:,:,itc) = tn(:,:,:,itc)
  where(ti(:,:,:,itc) < tn(:,:,:,itc)) ti(:,:,:,itc) = tn(:,:,:,itc)

end subroutine

!> Aborts the run if any state-vector field contains NaN values after the analysis step.
subroutine no_nan_constraint

  ! intern
  use state_module, only: state_vector, levX0, levX1, latX0, latX1, lonX0, lonX1

  ! tie-gcm
  use fields_module, only: f4d, itc, levd0, levd1, lond0, lond1, latd0, latd1

  implicit none

  integer :: i
  integer :: n_nans
  logical, dimension( levd0:levd1, lond0:lond1, latd0:latd1 ):: mask
  logical :: contains_nan

  contains_nan = .false.

  do i=1,size(state_vector%fd_idx,dim=1)
    mask = .false.
    where(isnan(f4d(state_vector%fd_idx(i))%data(:,:,:,itc))) mask = .true.
    n_nans= count(mask(levX0:levX1,lonX0:lonX1, latX0:latX1))
    if(n_nans > 0) then
      contains_nan = .true.
      write(*,"(3a,i6,a)") "ERROR. ", f4d(state_vector%fd_idx(i))%short_name, &
                           " contains ", n_nans, "  NAN after analysis step"
    end if
  end do

  if(contains_nan .eqv. .true.)then
    call shutdown("state vector contains NANs")
  end if

end subroutine

!> Clamps atomic oxygen ion (O+) density to a physically plausible range.
subroutine atomic_oxygen_ion_density_constraint

  ! tie-gcm
  use fields_module, only: itc, op, op_nm

  implicit none

  if(cfg_log%verbose_level>0) write(*,*) 'constraint atomic oxygen ion density'

  ! max: 1.2 * max value in TIE-GCM hallooween storm run rounded to two significant figures
  call limit_values(op(:,:,:,itc),1E-12, 4.0E+6)
  call limit_values(op_nm(:,:,:,itc),1E-12, 4.0E+6)

end subroutine

!> Clamps molecular oxygen ion (O2+) density to a physically plausible range.
subroutine molecular_oxygen_ion_density_constraint

  ! tie-gcm
  use fields_module, only: itc, o2p

  implicit none

  if(cfg_log%verbose_level>0) write(*,*) 'constraint molecular oxygen ion density'

  ! max: 1.2 * max value in TIE-GCM hallooween storm run rounded to two significant figures
  call limit_values(o2p(:,:,:,itc),1E-12, 2.9E+5)


end subroutine

!> Rescales ion densities (or overrides Ne) so the ionosphere is electrically neutral, i.e. Ne equals the summed ion density.
subroutine electrically_neutral_constraint(ions_are_ref)

  ! tie-gcm
  use fields_module,&
      only: itc, levd0, levd1, lond0, lond1, latd0, latd1
  use fields_module, only: itc, op, o2p, nop, nplus, n2p, ne, op_nm
  use mpi_module,&
    only: mytid
  use params_module,&
    only: nlat, nlon, nlev

  ! intern
  use mod_parallel_pdaf, only: task_id

  ! extern
!   use array_print_module, only: printMat

  implicit none

  ! local
  real, dimension(levd0:levd1, lond0:lond1, latd0:latd1) :: C, Cint, s
  logical, dimension(levd0:levd1, lond0:lond1, latd0:latd1) :: mask
  integer :: number_failed_cells

  ! arguments
  logical, intent(in) :: ions_are_ref

  C = op(:,:,:,itc)+o2p(:,:,:,itc)+nop(:,:,:)+nplus(:,:,:)+n2p(:,:,:) ! units: cm-3


  ! NE is given on interfaces
  ! ion densities are located on midpoints

  ! to interfaces
  Cint(levd0,:,:) = 1.5*C(levd0,:,:)-0.5*C(levd0+1,:,:)
  Cint(levd0+1:levd1,:,:) = 0.5*(  Cint(levd0+1:levd1,:,:)+ C(levd0:levd1-1,:,:) )

  mask = .false.
  if(ions_are_ref) then
    where( abs(Cint-ne(:,:,:,itc)) > .1 )
      mask = .true.
      ne(:,:,:,itc) = Cint
    end where
  else
    s = ne(:,:,:,itc)/Cint
    where( abs(Cint-ne(:,:,:,itc)) > .1 )
      mask = .true.

      op(:,:,:,itc) = op(:,:,:,itc)*s
      o2p(:,:,:,itc) = o2p(:,:,:,itc)*s
      nop(:,:,:) = nop(:,:,:)*s
      nplus(:,:,:) = nplus(:,:,:)*s
      n2p(:,:,:) = n2p(:,:,:)*s

    end where

    ! also compute for op_nm
    C = op_nm(:,:,:,itc)+o2p(:,:,:,itc)+nop(:,:,:)+nplus(:,:,:)+n2p(:,:,:)
    Cint(levd0,:,:) = 1.5*C(levd0,:,:)-0.5*C(levd0+1,:,:)
    Cint(levd0+1:levd1,:,:) = 0.5*(  Cint(levd0+1:levd1,:,:)+ C(levd0:levd1-1,:,:) )

    where( abs(Cint-ne(:,:,:,itc)) > .1 )
      op_nm(:,:,:,itc) = op_nm(:,:,:,itc) * ne(:,:,:,itc)/Cint
    end where
  end if

  call count_failed_cells(mask, number_failed_cells)

  if( mytid .eq. 0 ) then
    if (number_failed_cells .gt. 0) then
    write(*,*) 'ensemble member ', task_id, ' WARNING  ',&
                number_failed_cells, '/', nlat*nlon*nlev ,' cells', &
                ' are not electrically neutral'
    end if

  end if


end subroutine

!> Counts true entries of mask on the local subdomain and MPI-reduces the sum to rank 0.
subroutine count_failed_cells(mask, number_failed_cells)

  ! extern
  use mpi_f08

  ! tie-gcm
  use fields_module,&
      only:  levd0, levd1, lond0, lond1, latd0, latd1
  use mpi_module,&
      only: lon0, lon1, lat0, lat1, TIEGCM_WORLD, handle_mpi_err, mytid

  implicit none

  ! arguments
  logical, dimension(levd0:levd1, lond0:lond1, latd0:latd1), intent(in) :: mask
  integer, intent(out) :: number_failed_cells

  ! local
  integer :: reduced
  integer :: ierr

  number_failed_cells = count(mask(:,lon0:lon1, lat0:lat1))

  call MPI_Reduce(number_failed_cells, &
                  reduced, &
                  1, &
                  MPI_INTEGER, MPI_SUM, &
                  0, TIEGCM_WORLD, &
                  ierr)
  if (ierr /= MPI_SUCCESS) call  handle_mpi_err(ierr,'MPI reduce number of constrained cells')

  if( mytid .eq. 0 ) then
    number_failed_cells = reduced
  end if
end subroutine

!> Limits the relative change in the temporal gradient (rate of change between the previous and current time step) of each field.
subroutine gradient_constraint(fields,relative_limit)

  ! extern

  ! tie-gcm
  use fields_module,&
    only: f4d, itp, itc
  use fields_module, only: levd0, levd1, lond0, lond1, latd0, latd1

  ! intern
  use character_routines_module, only: string_ends_with

  implicit none

  ! arguments
  character(len=16), dimension(:), intent(in) :: fields
  real, intent(in) :: relative_limit

  ! local
  real, dimension( levd0:levd1, lond0:lond1, latd0:latd1 ):: grad_curr
  real, dimension( levd0:levd1, lond0:lond1, latd0:latd1 ):: grad_prev

  integer :: idx, idx_nm
  integer :: i

  do i=1, size(fields)
    select case(trim(fields(i)))
      case ("NE","O2P","N2D","TE","TI","OMEGA","POTEN")
        cycle
      case default
        if(string_ends_with(trim(fields(i)),'_NM').eqv..false.) then

          if(cfg_log%verbose_level>0) write(*,*) 'applying gradient constraint to', trim(fields(i)), &
                     '. Max relative change in gradient is: ', relative_limit

          idx = findloc(f4d%short_name,trim(fields(i)), dim=1)
          if(idx==0)then
            write(*,*) 'ERROR could not find ', trim(fields(i))
          end if

          idx_nm = findloc(f4d%short_name,trim(fields(i))//'_NM', dim=1)

          if(idx_nm==0)then
            write(*,*) 'ERROR could not find ', trim(fields(i))//'_NM'
          end if

          grad_curr = calc_gradient_f4d(f4d(idx)%data, f4d(idx_nm)%data, itc)
          grad_prev = calc_gradient_f4d(f4d(idx)%data, f4d(idx_nm)%data, itp)

          where(grad_curr<(1/relative_limit)*grad_prev)
            f4d(idx_nm)%data(:,:,:,itc) = f4d(idx)%data(:,:,:,itc) - (1/relative_limit)*grad_prev
          elsewhere(grad_curr>relative_limit*grad_prev)
            f4d(idx_nm)%data(:,:,:,itc) = f4d(idx)%data(:,:,:,itc) - relative_limit*grad_prev
          end where
        end if
    end select
  end do

end subroutine

!> Computes the temporal gradient of a TIE-GCM field, by finite-differencing its current-step data against its _NM (previous-timestep) counterpart at time index itx.
function calc_gradient_f4d(field, field_nm, itx) result(grad)

  use quantity_info_module, only: calc_gradient

  use fields_module, only: levd0, levd1, lond0, lond1, latd0, latd1

  implicit none

  ! arguments
  real, dimension(levd0:levd1, lond0:lond1, latd0:latd1, 2), intent(in) :: field
  real, dimension(levd0:levd1, lond0:lond1, latd0:latd1, 2), intent(in) :: field_nm
  integer, intent(in) :: itx

  ! result
  real, dimension(levd0:levd1, lond0:lond1, latd0:latd1) :: grad

  grad = calc_gradient(field(:,:,:,itx), field_nm(:,:,:,itx))

end function

! subroutine neutral_wind_constraint
!
!   ! tie-gcm
!   use fields_module, only: un, vn, un_nm, vn_nm, itc
!
!   call limit_values(un(:,:,:,itx),-1E+6,1E+6)
!
! end subroutine

! subroutine limit_update_by_previous(max_factor)
!
!   ! intern
!   use state_module, only: state_vector
!
!   ! tie-gcm
!   use fields_module, only: f4d, itc, itp
!
!   implicit none
!
!   real, intent(in) :: max_factor
!
!   ! local
!   integer :: i
!   real, dimension(:,:,:), pointer :: previous
!   real, dimension(:,:,:), pointer :: current
!
!   do i=1,size(state_vector%fd_idx,dim=1)
!     current => f4d(state_vector%fd_idx(i))%data(:,:,:,itc))
!     previous => f4d(state_vector%fd_idx(i))%data(:,:,:,itp))
!   end do
!
! end subroutine

!> Clamps array A element-wise to the given optional lower and/or upper bound.
subroutine limit_values(A, lower_limit, upper_limit)
  implicit none

  real, dimension(:,:,:), intent(inout) :: A
  real, optional, intent(in) :: lower_limit
  real, optional, intent(in) :: upper_limit

  if(present(lower_limit))  where( A < lower_limit ) A = lower_limit
  if(present(upper_limit))  where( A > upper_limit ) A = upper_limit

end subroutine


end module constraints
