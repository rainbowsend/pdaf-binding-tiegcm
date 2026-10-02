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
! Computes ensemble statistical moments (mean, variance, skewness, kurtosis) via MPI reduction over the ensemble.

module mpi_moments_module

implicit none

! this module provides subroutines for calculating statistical moments from data within COMM_Couple
!
! PDAF_diag_compute_moments cannot be used here, since it expects the whole
! ensemble to be available locally, i.e. ens(dim_p,dim_ens). Here every model
! task holds a single member instead, hence the reduction runs over
! COMM_couple, which decomposes the ensemble. The MPI routines of PDAF_diag
! reduce over COMM_filter, which decomposes the state, and thus solve the
! opposite problem. In addition the quantities are diagnostic ones such as ZG
! or VTEC, which are not part of the state vector and therefore unknown to
! PDAF.
!
! ATTENTION PDAF_unbiased_moments_from_summed_residuals below is a copy of the
! routine of the same name in PDAF_diag. PDAF declares it PRIVATE, so it
! cannot be imported. Keep both in sync.
!
! TODO ask upstream to make PDAF_unbiased_moments_from_summed_residuals and
! PDAF_biased_moments_from_summed_residuals public in PDAF_diag, so that the
! copy below can be dropped. Both routines originate from this project and
! were contributed to PDAF, only the declarations in PDAF_diag.F90 have to be
! removed:
!   PRIVATE PDAF_unbiased_moments_from_summed_residuals
!   PRIVATE PDAF_biased_moments_from_summed_residuals
! They are useful on their own whenever the summed exponentiated residuals are
! obtained differently than by PDAF_diag_compute_moments, here by an MPI
! reduction over the ensemble.

! allocates an array with an additional dimension for storing statistical moments
! (...,1) mean
! (...,2) variance
! (...,3) skewness
! (...,4) excess kurtosis
interface allocate_and_calc_moments_mpi
    procedure allocate_and_calc_moments_mpi_scalar,&
              allocate_and_calc_moments_mpi_1d, &
              allocate_and_calc_moments_mpi_2d, &
              allocate_and_calc_moments_mpi_3d
end interface allocate_and_calc_moments_mpi

! only calculates the mean
interface compute_mean_mpi
    procedure compute_mean_mpi_scalar,&
              compute_mean_mpi_mat_1d,&
              compute_mean_mpi_mat_2d,&
              compute_mean_mpi_mat_3d
end interface compute_mean_mpi

contains

 !> Allocates the moments output array on the filter root and computes the
 !! MPI-reduced statistical moments (up to order kmax) of a scalar across model tasks.
 subroutine allocate_and_calc_moments_mpi_scalar(x,kmax,moments)

    use mod_parallel_pdaf, only: filterpe

    implicit none
    real, intent(in) :: x
    integer, intent(in) :: kmax
    real, dimension(:), allocatable, intent(out) :: moments ! only on root (kmax)

    ! local
    real, dimension(1) :: x_
    real, dimension(:,:), allocatable :: moments_

    if(filterpe) then
      allocate (moments(kmax))
      moments = 0.0
      allocate(moments_(1,kmax))
    else
      allocate (moments(0))
      allocate (moments_(0,0))
    end if

    x_(1) = x
    call calc_moments_mpi(x_,kmax,moments_)
    moments = moments_(1,:)

    deallocate(moments_)

  end subroutine

  !> Allocates the moments output array on the filter root and computes the
  !! MPI-reduced statistical moments (up to order kmax) of a 1D array across model tasks.
  subroutine allocate_and_calc_moments_mpi_1d(x,kmax,moments)

    use mod_parallel_pdaf, only: filterpe

    implicit none
    real, dimension(:), intent(in) :: x
    integer, intent(in) :: kmax
    real, dimension(:,:), allocatable, intent(out) :: moments ! only on root (:,kmax)

    if(filterpe) then
      allocate (moments(size(x,dim=1), kmax))
      moments = 0.0
    else
      allocate (moments(0,0))
    end if
    call calc_moments_mpi(x,kmax,moments)

  end subroutine

  !> Allocates the moments output array on the filter root and computes the
  !! MPI-reduced statistical moments (up to order kmax) of a 2D array (flattened
  !! before reduction) across model tasks.
  subroutine allocate_and_calc_moments_mpi_2d(x,kmax,moments)

    use array_mapping_module, only: flatten
    use mod_parallel_pdaf, only: filterpe

    implicit none
    real, dimension(:,:), intent(in), contiguous, target :: x
    integer, intent(in) :: kmax
    real, dimension(:,:,:), allocatable, intent(out), target :: moments ! only on root (:,:,kmax)

    real, dimension(:), pointer, contiguous :: x_flat
    real, dimension(:,:), pointer, contiguous :: moments_flat

    x_flat => flatten(x)

    if(filterpe) then
      allocate (moments(size(x,dim=1), size(x,dim=2), kmax))
      moments = 0.0
      moments_flat(1:size(x),1:kmax) => moments(:,:,:)
    else
      allocate (moments(0,0,0))
      moments_flat(1:1,1:1) => moments
    end if
    call calc_moments_mpi(x_flat,kmax,moments_flat)

  end subroutine

  !> Allocates the moments output array on the filter root and computes the
  !! MPI-reduced statistical moments (up to order kmax) of a 3D array (flattened
  !! before reduction) across model tasks.
  subroutine allocate_and_calc_moments_mpi_3d(x,kmax,moments)

    use array_mapping_module, only: flatten
    use mod_parallel_pdaf, only: filterpe

    implicit none

    ! arguments
    real, dimension(:,:,:), intent(in), contiguous, target :: x
    integer, intent(in) :: kmax
    real, dimension(:,:,:,:), allocatable, intent(out), target :: moments ! only on root (:,:,:,kmax)

    ! local
    real, dimension(:), pointer, contiguous :: x_flat
    real, dimension(:,:), pointer, contiguous :: moments_flat

    x_flat => flatten(x)

    if(filterpe) then
      allocate (moments(size(x,dim=1), size(x,dim=2), size(x,dim=3), kmax))
      moments = 0.0
      moments_flat(1:size(x),1:kmax) => moments(:,:,:,:)
    else
      allocate (moments(0,0,0,0))
      moments_flat(1:1,1:1) => moments
    end if

    call calc_moments_mpi(x_flat,kmax,moments_flat)

  end subroutine

  !> Computes up to the kmax-th central moment (mean, variance, skewness, excess
  !! kurtosis) of x by MPI-reducing over the model-task communicator; results are
  !! only valid on filter PEs.
  subroutine calc_moments_mpi(x,kmax,moments)

    ! extern
    use mpi_f08

    ! intern
    use mod_parallel_pdaf, only: COMM_couple, filterpe
    use mod_parallel_pdaf, only: n => n_modeltasks
!     use array_print_module, only: printMat

    implicit none

    ! arguments
    real, dimension(:), intent(in) :: x
    integer, intent(in) :: kmax
    real, dimension(:,:), intent(out) :: moments ! need to be allocated only on root (x,kmax)

    ! local
    real, dimension(:), allocatable :: mean
    real, dimension(:,:), allocatable :: residuals ! [r r**2 r**3 ... r**kmax]

    integer :: kmax_

    integer :: i
    integer :: ierr
    integer :: dim_x

    if(kmax>n) then
      write(*,'(a,i1,a)') 'WARNING not enough samples to compute ', kmax, '-th moment'
      kmax_ = n
    else
      kmax_ = kmax
    end if

    dim_x = size(x)

    allocate(mean(dim_x))

    ! first moment (mean)
    call MPI_Reduce( x, mean, dim_x,   &
                    MPI_REAL8, MPI_SUM,  0,    &
                    COMM_couple, ierr);

    mean = mean/n
    if(filterpe) then
      moments(:,1) = mean
    end if

    if( kmax_ > 1 ) then
      call MPI_Bcast( mean, dim_x,   &
                        MPI_REAL8,  0,    &
                        COMM_couple,  ierr);

      ! compute residuals and power of residuals [r**1 r**2 r**3 ... r**kmax]
      ! overflow in kurtosis calculation occures if residulal > 10**77 (assuming E+308 is largest number)
      allocate(residuals(dim_x,1:kmax_))

      residuals(:,1) = x - mean
      do i = 2, kmax_
        residuals(:,i) = residuals(:,i-1)*residuals(:,1)
      end do

      call MPI_Reduce( residuals(:,2:kmax_), moments(:,2:kmax_), dim_x*(kmax_-1),   &
                      MPI_REAL8, MPI_SUM,  0,    &
                      COMM_couple,  ierr);

      if(filterpe) then
        call PDAF_unbiased_moments_from_summed_residuals(n,&
                                  dim_x,&
                                  kmax_,&
                                  moments(:,1:kmax_),&
                                  moments(:,1:kmax_))
      end if

      deallocate(residuals)
    end if

    deallocate(mean)


  end subroutine

 !> Computes the MPI-reduced mean of a scalar across model tasks.
 subroutine compute_mean_mpi_scalar(x,mean)

    ! extern
    use mpi_f08

    ! intern
    use mod_parallel_pdaf, only: COMM_couple,  n_modeltasks

    implicit none

    ! arguments
    real, intent(in) :: x
    real, intent(out) :: mean

    ! local
    integer :: ierr

    if(n_modeltasks == 1) then
      mean = x
    else

      call MPI_reduce( x, mean, 1,   &
                  MPI_REAL8, MPI_SUM,  0,    &
                  COMM_couple,  ierr);

      mean = mean/n_modeltasks;

    end if

  end subroutine compute_mean_mpi_scalar

  !> Computes the elementwise MPI-reduced mean of a 1D array across model tasks.
  subroutine compute_mean_mpi_mat_1d(x,mean)

    ! extern
    use mpi_f08

    ! intern
    use mod_parallel_pdaf, only: COMM_couple,  n_modeltasks

    implicit none

    ! arguments
    real, dimension(:), intent(in) :: x
    real, dimension(:), intent(out) :: mean

    ! local
    integer :: ierr

    if(n_modeltasks == 1) then
      mean = x
    else

      call MPI_reduce( x, mean, size(x),   &
                  MPI_REAL8, MPI_SUM,  0,    &
                  COMM_couple,  ierr);

      mean = mean/n_modeltasks;

    end if

  end subroutine compute_mean_mpi_mat_1d

  !> Computes the elementwise MPI-reduced mean of a 2D array across model tasks.
  subroutine compute_mean_mpi_mat_2d(x,mean)

    ! extern
    use mpi_f08

    ! intern
    use mod_parallel_pdaf, only: COMM_couple,  n_modeltasks

    implicit none

    ! arguments
    real, dimension(:,:), intent(in) :: x
    real, dimension(:,:), intent(out) :: mean

    ! local
    integer :: ierr

    if(n_modeltasks == 1) then
      mean = x
    else

      call MPI_reduce( x, mean, size(x),   &
                  MPI_REAL8, MPI_SUM,  0,    &
                  COMM_couple,  ierr);

      mean = mean/n_modeltasks;

    end if

  end subroutine compute_mean_mpi_mat_2d

  !> Computes the elementwise MPI-reduced mean of a 3D array across model tasks.
  subroutine compute_mean_mpi_mat_3d(x,mean)

    ! extern
    use mpi_f08

    ! intern
    use mod_parallel_pdaf, only: COMM_couple,  n_modeltasks

    implicit none

    ! arguments
    real, dimension(:,:,:), intent(in) :: x
    real, dimension(:,:,:), intent(out) :: mean

    ! local
    integer :: ierr

    if(n_modeltasks == 1) then
      mean = x
    else

      call MPI_reduce( x, mean, size(x),   &
                  MPI_REAL8, MPI_SUM,  0,    &
                  COMM_couple,  ierr);

      mean = mean/n_modeltasks;

    end if

  end subroutine compute_mean_mpi_mat_3d


!--------------------------------------------------------------------------
!> Computes the unbiased estimator for mean, variance, skewness, and excess
!! kurtosis from the sum of exponentiated residulals
!!
!! Computes the unbiased estimator for mean, variance, skewness, and excess
!! kurtosis from the sum of exponentiated residulals. You can perform an inplace
!! moment calculation by using the same input for sum_expo_resid and moments
!!
!! __Revision history:__
!! * 2023-08 - Armin Corbin - original code for tiegcm-pdaf
!! * 2025-03 - Armin Corbin - ported for PDAF 3
!!
!! TODO copy of the PRIVATE routine of the same name in PDAF_diag, see the
!! header of this module. Delete it once it is public in PDAF.
SUBROUTINE PDAF_unbiased_moments_from_summed_residuals(dim_ens, dim_p, kmax, sum_expo_resid, moments)

  IMPLICIT NONE

  ! *** Arguments ***
  INTEGER, INTENT(IN) :: dim_ens  !< number of ensemble members/samples
  INTEGER, INTENT(IN) :: dim_p    !< local size of the state
  INTEGER, INTENT(IN) :: kmax     !< maximum order of central moment that is computed, maximum is 4
  REAL, INTENT(IN) :: sum_expo_resid(dim_p, kmax) ! sum of exponentiated residulals, first col is ignored
                                  !< [1st, sum(r**2), sum(r**3), ... ,sum(r**kmax)]
  REAL, INTENT(INOUT) :: moments(dim_p, kmax)     !  The columns contain the moments of the ensemble
                                  !< (mean, variance, skewness, excess kurtosis)

  ! unbiased estimator of variance
  ! k2 = sum(r**2)/(n-1)
  moments(:,2) = sum_expo_resid(:,2)/REAL(dim_ens-1)

  ! unbiased skewness
  IF(kmax>2) THEN
    WHERE(moments(:,2)/=0)
      ! G1 = k3/k2**(3/2)
      moments(:,3) = REAL(dim_ens)/REAL((dim_ens-1)*(dim_ens-2)) &
                    * sum_expo_resid(:,3)/moments(:,2)**(3./2.)
    ELSEWHERE
      moments(:,3) = 0.
    ENDWHERE
  END IF

  ! unbiased excess kurtosis
  IF(kmax>3) THEN
    ! G2 = k4/k2**2
    WHERE(moments(:,2)/=0)
      moments(:,4) = REAL(dim_ens*(dim_ens+1))/REAL((dim_ens-1)*(dim_ens-2)*(dim_ens-3)) &
                      * sum_expo_resid(:,4)/moments(:,2)**2
      moments(:,4) = moments(:,4) - 3.*REAL((dim_ens-1)**2)/REAL(((dim_ens-2)*(dim_ens-3)))
    ELSEWHERE
      moments(:,4) = 0.
    ENDWHERE
  END IF

END SUBROUTINE PDAF_unbiased_moments_from_summed_residuals

end module
