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
! Initializes/manages the ensemble: per-member source states, cross-rank gathering, pseudo-inverse-based calibration updates.

module ensemble_module
! ensemble initalization

  ! intern
  use array_print_module, only: printMat
  use state_module,  only: state_vector

  implicit none

  logical :: init_ens_from_first_state

  real, dimension(:,:), allocatable, save :: state_0 ! first state of all ensemble (local_state_size, ensemble size)

contains

!> Deallocates the module-level state_0 array holding the initial ensemble states.
SUBROUTINE deallocate_ensemble
    implicit none
     if (allocated (state_0)) deallocate (state_0)
END SUBROUTINE deallocate_ensemble

!> Builds this ensemble task's per-member source-file path from the configured ensemble source path/name and task_id, and aborts if that file doesn't exist.
subroutine set_source()

  ! tie-gcm
  use input_module, only: source

  ! intern
  use configuration, only: cfg_ensemble
  use mod_parallel_pdaf, only: task_id

  implicit none

  logical :: is_existent 

  source = ''
  write(source,'(2a,i0.4,2a)') trim(cfg_ensemble%source_path), '/ens_', task_id, '_', trim(cfg_ensemble%source_name)

  write(*,*) '! Ignoring source in input file, replacing it wih: ', source

  inquire( file=trim(source), exist=is_existent ) 
  if ( is_existent .neqv. .true. ) then
    call shutdown('source file: ' // trim(source) // 'does not exist')
  end if

end subroutine

!> Fills state_0 (the ensemble's initial/first state, gathered onto the filter PE) from
!! each ensemble member's local state, initialized from that member's own primary history
!! file (see set_source) and optionally overlaid with calibration-parameter values. Must
!! run before PDAF_init, since PDAF's init_ens callback later just copies from state_0.
subroutine fill_state_0 ()

  ! extern
  use mpi_f08

  ! tie-gcm
  use fields_module, only: itp ! sources are read into itp
  use mpi_module, only: mytid, handle_mpi_err

  ! intern
  use configuration, only: cfg_calibration
  use mod_parallel_pdaf, only: filterpe, n_modeltasks, COMM_couple

  implicit none

  real, allocatable, dimension(:) :: local_state
  integer :: dim_state_p
  integer :: MPIerr

  dim_state_p = state_vector%map%size_R(mytid)

  if( filterpe )then
    allocate( state_0(dim_state_p,n_modeltasks) )
  end if

  allocate( local_state(dim_state_p) )

  call state_vector%fill_state_p(local_state, itp)

  if(cfg_calibration%apply) call read_calibration_ensemble(local_state)

  call MPI_Gather(local_state, dim_state_p, MPI_REAL8, &
                      state_0, dim_state_p, MPI_REAL8, &
                            0,                         &
                      COMM_couple, MPIerr)
  if (MPIerr /= MPI_SUCCESS) call handle_mpi_err(MPIerr,'fill_state_0')

  deallocate(local_state)

end subroutine

!> On the root PE of this model instance (mytid==0 within TIEGCM_WORLD, i.e. once per
!! ensemble member), overwrites the model inputs augmenting the state vector (being
!! co-estimated) with the corresponding values read from the perturbations file. Used
!! for initializing the augmented part of the state vector.
subroutine read_calibration_ensemble(local_state)

  ! tie-gcm
  use mpi_module, only: mytid

  ! intern
  use array_print_module,&
    only: printMat
  use model_parameter_handling_module,&
    only: model_parameters

  implicit none

  ! arguments
  real, dimension(:), target, intent(inout) :: local_state

  ! local
  integer :: i
  integer :: idx
  real, dimension(:), contiguous, pointer :: sub_state

  if(mytid==0)then
    do i= state_vector%idx_cal_0, state_vector%idx_cal_1
      idx = findloc( model_parameters%name,state_vector%map%fd_name(i),dim=1 )
      sub_state => local_state(state_vector%map%idx_R(i, 0)%begin_p : &
                               state_vector%map%idx_R(i, 0)%back_p)
      sub_state = model_parameters(idx)%val
!       write(*,*) 'initalizing parameter ',  model_parameters(idx)%name , ' as ', sub_state
    end do
  end if

end subroutine

!> Gathers the ensemble matrix ens_p, distributed row-wise (by subdomain) over the
!! filter communicator COMM_filter, into the full ensemble matrix ens on rank 0 of
!! COMM_filter (= rank_world 0). Uses custom interleaved MPI row datatypes.
subroutine gather_ensemble(dim_p, dim_ens, ens_p, ens)

  ! extern
  use mpi_f08

  ! tie-gcm
  use mpi_module, only: handle_mpi_err

  ! intern
  use mod_parallel_pdaf, only: COMM_filter, rank_world, local_ntask

  implicit none

  ! args
  integer, intent(in) :: dim_p, dim_ens
  real, dimension(dim_p, dim_ens), intent(in) :: ens_p
  real, dimension(:,:), allocatable, intent(out) :: ens

  ! local
  integer :: rank
  integer :: ierror
  integer :: dim_full_state

  integer, dimension(local_ntask) :: counts, displs

  TYPE(MPI_Datatype) :: row_p_tmp, row_tmp
  TYPE(MPI_Datatype) :: row_p, row

  integer :: i

  INTEGER(KIND=MPI_ADDRESS_KIND) :: lb, extent
  real :: r

  call MPI_Comm_rank(COMM_filter, rank, ierror)

  call MPI_Allgather(dim_p, 1, MPI_Integer, &
                     counts, 1, MPI_Integer, &
                     COMM_filter,ierror)

  dim_full_state = sum(counts)

  displs(1) = 0
  do i=2,local_ntask
    displs(i)=displs(i-1)+counts(i-1)
  end do

  if(rank_world==0)then
    if(allocated(ens)) deallocate(ens)
    allocate(ens(dim_full_state,dim_ens))
  else
    allocate(ens(0,0))
  end if

  ! Custom row types are needed because ens_p/ens are column-major with different
  ! leading dimensions (dim_p vs dim_full_state): a plain Gatherv would concatenate
  ! each rank's contiguous block byte-for-byte instead of interleaving rows by column.

  ! TODO save type and do not create it each call
  ! subtype for a single row of matrix
  CALL MPI_Type_vector(dim_ens, 1, dim_full_state, MPI_REAL8, row_tmp)
  CALL MPI_Type_vector(dim_ens, 1, dim_p, MPI_REAL8, row_p_tmp)

  ! Resize row type to make sure it is interleaved when repeated
  lb=0
  extent=sizeof(r)
  call MPI_Type_create_resized (row_tmp, lb, extent, row, ierror)
  call MPI_Type_create_resized (row_p_tmp, lb, extent, row_p, ierror)

  call MPI_Type_commit(row, ierror)
  call MPI_Type_commit(row_p, ierror)

  call MPI_Gatherv(ens_p, dim_p, row_p, &
                   ens, counts, displs, row, &
                   0, COMM_filter,ierror)
  if (ierror /= MPI_SUCCESS) call handle_mpi_err(ierror,'gather_ensemble')


  call MPI_Type_free(row,ierror);
  call MPI_Type_free(row_p,ierror);

  if(rank_world/=0)then
    deallocate(ens)
  end if


end subroutine

!> Gathers the full ensemble and computes its Moore-Penrose pseudo-inverse via SVD on rank_world 0, dropping singular values below a computed tolerance.
subroutine pseudo_inverse_ensemble_like(dim_p, dim_ens, ens_p, pinv)
  ! compute the pseude inverse using singular value decomposition
  !
  ! A = U S V'
  !
  ! INV = V S-1 U*
  !
  ! singular values smaller than a computed tolerance are ignored

  ! intern
  use default_lapack_interface, only: default_dgemm
  use mod_parallel_pdaf, only:  rank_world
  use netcdf_functionality, only: write_mat_2d_to_netcdf

  implicit none

  external :: dgesvd, dgemm, dscal

  !arguments
  integer, intent(in) :: dim_p, dim_ens
  real, dimension(dim_p, dim_ens), intent(in) :: ens_p
  real, dimension(:,:), allocatable, intent(out) :: pinv

  ! local
  real, dimension(:,:), allocatable :: ens

  integer :: M,N,LWORK
  integer :: INFO
  real, dimension(dim_ens) :: SVA
  real, dimension(:,:), allocatable :: U,VT
  real, dimension(:), allocatable :: WORK
  real, dimension(1) :: wq

  integer :: i

  integer :: significant_sv

  real :: tol

  ! TODO use scalapack for parallel computation of SVD

  call gather_ensemble(dim_p, dim_ens, ens_p, ens)

  if(rank_world==0)then

!     call write_mat_2d_to_netcdf(ens,'svd_debug_A.nc')

    M=size(ens,dim=1)
    N=dim_ens

    allocate(U(M,N))
    allocate(VT(N,N))

    write(*,*) 'computing SVD'

    ! first call querries optimal workspace size
    ! second call computes svd
    call dgesvd('S','A',M,N,ens,M,SVA,U,M,VT,N,wq,-1,info)

    LWORK=int(wq(1))
    allocate(WORK(LWORK))

    call dgesvd('S','A',M,N,ens,M,SVA,U,M,VT,N,WORK,LWORK,info)
    if(INFO/=0) write(*,*) "LAPACK Error. SVD computation in pseudo_inverse_ensemble_like"

    write(*,*) 'finished computing SVD'

!     call write_mat_2d_to_netcdf(U,'svd_debug_U.nc')
!     call write_mat_2d_to_netcdf(VT,'svd_debug_VT.nc')

    tol=M*epsilon(maxval(SVA))*1000

    ! invert singular values, if they are larger than the tolerance
    where (SVA<tol)
      SVA=0
    elsewhere
      SVA=1./SVA
    end where

    ! multiply each column of U with corresponding singular value
    ! here we calculate (S^-1 * U')'
    do i=1,N
      if(SVA(i)==0)then
        exit
      end if
      significant_sv=i
      call dscal(M,SVA(i),U(:,i),1)
    end do

    write(*,*) 'tolerance is ', tol,  ' -> number of significant singular values:', significant_sv, ' / ', dim_ens

    if(allocated(pinv)) deallocate(pinv)
    allocate(pinv(N,M))

    ! compute V * [S^-1 * U']
    ! We only need to compute the matrix multiplication for the first n rows of [S^-1 * U']
    ! and n cols of V corresponding to the non zero singular values
    !
    call default_dgemm(VT,U,pinv,tA=.true.,tB=.true.,K=significant_sv)

    ! TODO running the same compiled code twice yields slightly different results for SVD and pinv, probably by reason of some random generated values in svd calculation

!     call write_mat_2d_to_netcdf(pinv,'svd_debug_pinv.nc')

    deallocate(U)
    deallocate(VT)
    deallocate(WORK)

  end if

  if (allocated(ens)) deallocate(ens)

end subroutine

!> Updates the calibration-parameter block of the augmented ensemble via linear regression on the state-vector analysis increment (Bocquet et al. 2021, Eq. 21), scaled by a configured tapering factor.
subroutine update_calibration_parameters_lin_reg(dim_state_p, dim_cal, dim_ens, ens_analysis_p, ens_forecast_p)
!
! Input arguments are the augumented ensemble at the forecast and analysis step
!
! Augumented state vector: z = [x p]
!
! ATTENTION the calibration parameters p are appended to the state vector of the root.
! The other ranks do not hold any calibration parameters
!
! filter rank    0    1    2   ...  n-1
!               |x|  |x|   |x|      |x|
!               |x|  |x|   |x|      |x|
!               |x|  |x|   |x|      |x|
!               |x|  |x|   |x|      |x|
!               |p|
!
!

  ! intern
  use configuration, only: cfg_calibration
  use default_lapack_interface, only: default_dgemm
  use mod_parallel_pdaf, only: rank_filter
  use netcdf_functionality, only: write_mat_2d_to_netcdf

  implicit none

  ! arguments
  integer, intent(in) :: dim_state_p ! number of states in augumented vector
  integer, intent(in) :: dim_cal     ! number of calibration parameters in augumented vector
  integer, intent(in) :: dim_ens

  real, dimension(dim_state_p+dim_cal, dim_ens), target, intent(inout) :: ens_analysis_p
  real, dimension(dim_state_p+dim_cal, dim_ens), target, intent(in) :: ens_forecast_p

  ! local
  real, dimension(:,:), pointer :: ens_state_analysis_p
  real, dimension(:,:), pointer :: ens_state_forecast_p
  real, dimension(:,:), pointer :: ens_cal_forecast
  real, dimension(:,:), pointer :: ens_cal_analysis

  real, dimension(dim_state_p, dim_ens) :: ens_pert_state_forecast_p
  real, dimension(dim_cal, dim_ens)   :: ens_pert_cal_forecast

  real, dimension(dim_state_p) :: ens_mean_state_forecast_p
  real, dimension(dim_cal) :: ens_mean_cal_forecast

  real, dimension(dim_state_p, dim_ens) :: ens_state_delta_p

  real, dimension(:,:), allocatable :: ens_state_delta
  real, dimension(:,:), allocatable :: ens_pert_state_forecast_inv
  real, dimension(:,:), allocatable :: solution
  real, dimension(:,:), allocatable :: cal_update

  real :: norm

  write(*,*) 'updating global parameters'

  ! extract 'state' and 'dynamics' from augumented vector

  ens_state_analysis_p => ens_analysis_p(1:dim_state_p,:)
  ens_state_forecast_p => ens_forecast_p(1:dim_state_p,:)

  if (rank_filter == 0) then
    ens_cal_forecast  => ens_forecast_p(dim_state_p+1:dim_state_p+dim_cal,:)
    ens_cal_analysis  => ens_analysis_p(dim_state_p+1:dim_state_p+dim_cal,:)
  else
    ! this causes compiler warning ‘ens_cal_forecast.dim[1].stride’ may be used uninitialized in this function
    ens_cal_forecast  => null()
    ens_cal_analysis  => null()
  end if

  norm = sqrt(float(dim_ens)-1.)

  ! compute normalized perturbations
  ens_mean_state_forecast_p = sum(ens_state_forecast_p,dim=2)/dim_ens
  ens_pert_state_forecast_p = (ens_state_forecast_p-spread(ens_mean_state_forecast_p,dim=2,ncopies=dim_ens))/norm

  if (rank_filter == 0) then
    ens_mean_cal_forecast = sum(ens_cal_forecast,dim=2)/dim_ens
    ens_pert_cal_forecast = (ens_cal_forecast-spread(ens_mean_cal_forecast,dim=2,ncopies=dim_ens))/norm
  end if


  ! ---------------------------------------------------------------
  ! Eq. 21 in Bocquet et al. (2021)
  ! https://www.aimsciences.org/article/doi/10.3934/fods.2020015
  ! ---------------------------------------------------------------

  ens_state_delta_p = ens_state_analysis_p - ens_state_forecast_p

  call gather_ensemble(dim_state_p, dim_ens, ens_state_delta_p, ens_state_delta)

!   if (rank_filter == 0) then
!     call write_mat_2d_to_netcdf(ens_state_delta,"debug_cal_delta")
!     call write_mat_2d_to_netcdf(ens_pert_cal_forecast,"debug_cal_perturb")
!   end if

  ! TODO dgelsd, dgelss dgels
  call pseudo_inverse_ensemble_like(dim_state_p,&
                                    dim_ens,&
                                    ens_pert_state_forecast_p,&
                                    ens_pert_state_forecast_inv)

  if (rank_filter == 0) then

    ! S = inv(X_fx) * (E_ax - E_fx)
    allocate(solution(dim_ens,dim_ens))
    call default_dgemm(ens_pert_state_forecast_inv,ens_state_delta,solution)

    ! compute update for calibration parameters
    ! update = X_fp * S
    allocate(cal_update(dim_cal,dim_ens))
    call default_dgemm(ens_pert_cal_forecast, &
                       solution,&
                       cal_update)


    ! E_ap = E_fp + tapering * update

    ens_cal_analysis = ens_cal_forecast  + cfg_calibration%localization_tapering * cal_update

!     call write_mat_2d_to_netcdf(cal_update,"debug_cal_update")

  end if

  if(allocated(ens_state_delta)) deallocate(ens_state_delta)
  if(allocated(ens_pert_state_forecast_inv)) deallocate(ens_pert_state_forecast_inv)
  if(allocated(solution)) deallocate(solution)
  if(allocated(cal_update)) deallocate(cal_update)

end subroutine


end module ensemble_module
