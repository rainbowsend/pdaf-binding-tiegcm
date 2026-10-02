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
! Fortran wrappers around LAPACK/BLAS routines (dgemm, dposv, dgesv, dgeev) used by ensemble/calibration code.

module default_lapack_interface

contains

!> Wrapper around LAPACK dgemm: computes C = alpha*op(A)*op(B) + beta*C. tA/tB default to
!! .false. (not transposed) if omitted; M/N/K default to the corresponding dimensions of
!! A/B (consistent with the resulting tA_/tB_), otherwise are taken as given.
subroutine default_dgemm(A,B,C,alpha,beta,tA,tB,M,N,K)

  implicit none

  external :: dgemm

  ! arguments
  real, dimension(:,:), intent(in) :: A
  real, dimension(:,:), intent(in) :: B
  real, dimension(:,:), intent(out) :: C

  real, intent(in), optional :: alpha
  real, intent(in), optional :: beta

  logical, intent(in), optional :: tA
  logical, intent(in), optional :: tB

  integer, intent(in), optional :: M
  integer, intent(in), optional :: N
  integer, intent(in), optional :: K

  ! local
  real, parameter :: s=1.0

  real :: alpha_, beta_
  logical :: tA_, tB_
  integer :: M_, N_, K_

  character :: transA, transB

  if(present(alpha))then
    alpha_ = alpha
  else
    alpha_ = s
  end if

  if(present(beta))then
    beta_ = beta
  else
    beta_ = s
  end if

  if(present(tA))then
    tA_ = tA
  else
    tA_ = .false.
  end if

  if(present(tB))then
    tB_ = tB
  else
    tB_ = .false.
  end if

  if(tA_) then
    transA = 'T'
    if(present(M).eqv..false.) M_=size(A,dim=2)
  else
    transA = 'N'
    if(present(M).eqv..false.) M_=size(A,dim=1)
  end if

  if(tB_) then
    transB = 'T'
    if(present(N).eqv..false.) N_=size(B,dim=1)
    if(present(K).eqv..false.) K_=size(B,dim=2)
  else
    transB = 'N'
    if(present(N).eqv..false.) N_=size(B,dim=2)
    if(present(K).eqv..false.) K_=size(B,dim=1)
  end if

  if(present(M)) M_ = M
  if(present(N)) N_ = N
  if(present(K)) K_ = K

  C=0

  call dgemm(transA,&
             transB, &
             M_, &
             N_, &
             K_, &
             alpha_, &
             A, &
             size(A,dim=1), &
             B, &
             size(B,dim=1), &
             beta_, &
             C, &
             size(C,dim=1))

end subroutine

!> Wrapper around LAPACK dposv: solves A*X = B for a symmetric positive-definite A via Cholesky factorization, returning X in B.
subroutine default_dposv(UPLO,A,B,n,nRHS)

  ! arguments
  real, dimension(:,:), intent(inout) :: A ! On exit, if INFO = 0, the factor U or L from the Cholesky factorization A = U**T*U or A = L*L**T.
  real, dimension(:,:), intent(inout) :: B ! On exit, if INFO = 0, the N-by-NRHS solution matrix X

  character(len=1), intent(in) :: UPLO

  integer, intent(in), optional :: n
  integer, intent(in), optional :: nRHS

  ! local
  integer :: n_
  integer :: nRHS_

  integer :: INFO

  if(present(n).eqv..false.) then
    n_=size(A,dim=1)
  else
    n_=n
  end if

  if(present(nRHS).eqv..false.) then
    nRHS_=size(B,dim=2)
  else
    nRHS_=nRHS
  end if

  call dposv(UPLO, &
             n_, &
             nRHS_, &
             A, &
             size(A,dim=1), &
             B, &
             size(B,dim=1), &
             INFO)
   IF (info /= 0) THEN
      WRITE (*, *) '!!! dposv error code ', info
  ENDIF

end subroutine

!> Wrapper around LAPACK dgesv: solves A*X = B for a general N-by-N matrix A via LU
!! factorization with partial pivoting, returning X in B. Overwrites A in place with its
!! LU factors.
subroutine default_dgesv(A,b,n,nRHS)

  ! arguments
  real, dimension(:,:), intent(inout) :: A
  real, dimension(:,:), intent(inout) :: b

  integer, intent(in), optional :: n
  integer, intent(in), optional :: nRHS

  ! local
  integer :: n_
  integer :: nRHS_

  integer, dimension(:), allocatable :: IPIV

  integer :: INFO

  if(present(n).eqv..false.) then
    n_=size(A,dim=1)
  else
    n_=n
  end if

  if(present(nRHS).eqv..false.) then
    nRHS_=size(b,dim=2)
  else
    nRHS_=nRHS
  end if

  allocate(IPIV(n_))

  call DGESV( n_, nRHS_, A, size(A,dim=1), IPIV, B, size(b,dim=1), INFO )

  deallocate(IPIV)
end subroutine

!> Solves A*x = b like default_dgesv, but leaves A and b unmodified by operating on internal copies and returning the solution in x.
subroutine default_dgesv_not_inplace(A,x,b)

  ! arguments
  real, dimension(:,:), intent(in) :: A
  real, dimension(:,:), intent(in) :: b
  real, dimension(:,:), intent(out) :: x

  real, dimension(:,:), allocatable :: A_
  real, dimension(:,:), allocatable :: b_

  allocate(A_,source=A)
  allocate(b_,source=b)
  call default_dgesv(A_,b_)
  x = b_
  deallocate(A_)
  deallocate(b_)

end subroutine

!> Wrapper around LAPACK dgeev: computes the eigenvalues (and, if requested, left/right eigenvectors) of a general real matrix A, auto-sizing the workspace query internally.
subroutine default_DGEEV(A,eigenvalues,n,Vleft,Vright)

  ! arguments
  real, dimension(:,:), intent(inout) :: A
  integer, intent(in), optional :: n
  complex, dimension(:), intent(out) :: eigenvalues
  real, dimension(:), intent(out), target, optional :: Vright
  real, dimension(:), intent(out), target, optional :: Vleft

  ! local
  character(len=1) :: jobvl
  character(len=1) :: jobvr

  integer :: n_
  real, dimension(:), allocatable :: Wreal
  real, dimension(:), allocatable :: Wimag
  real, dimension(:), allocatable :: Work

  real, dimension(:), pointer :: Vleft_
  real, dimension(:), pointer :: Vright_

  real :: optimal_size(1)

  integer :: lwork
  integer :: info
  integer :: ldvl
  integer :: ldvr

  integer :: i

  if(present(Vleft)) then
    jobvl = 'V'
    ldvl = size(Vleft,dim=1)
    Vleft_ => Vleft
  else
    jobvl = 'N'
    ldvl = 1
    Vleft_ => null()
  end if

  if(present(Vright)) then
    jobvr = 'V'
    ldvr = size(Vright,dim=1)
    Vright_ => Vright
  else
    jobvr = 'N'
    ldvr = 1
    Vright_ => null()
  end if

  if(present(n).eqv..false.) then
    n_=size(A,dim=2)
  else
    n_=n
  end if

  allocate(Wreal(n_))
  allocate(Wimag(n_))

  call dgeev(jobvl, &
             jobvr, &
              n_, &
              A, &
              size(A,dim=1), &
              Wreal, &
              Wimag, &
              Vleft_, &
              ldvl, &
              Vright_, &
              ldvr, &
              optimal_size, &
              -1, &
              info)

  lwork=int(optimal_size(1))

  allocate(work(lwork))
  call dgeev(jobvl, &
            jobvr, &
            n_, &
            A, &
            size(A,dim=1), &
            Wreal, &
            Wimag, &
            Vleft_, &
            ldvl, &
            Vright_, &
            ldvr, &
            work, &
            lwork, &
            info)

  do i=1,n_
    eigenvalues(i) = COMPLEX(Wreal(i),Wimag(i))
  end do


  deallocate(Wreal)
  deallocate(Wimag)
  deallocate(work)

end subroutine

!> Eigenvalues-only variant of default_DGEEV for a fixed n-by-n matrix that caches the LAPACK workspace size across calls with the same n.
subroutine default_DGEEV_nxn(A,n,eigenvalues)

  ! arguments
  real, dimension(n,n), intent(inout) :: A
  integer, intent(in) :: n
  complex, dimension(n), intent(out) :: eigenvalues


  ! local
  real, dimension(n) :: Wreal
  real, dimension(n) :: Wimag

  real, dimension(:), allocatable :: Work

  real, dimension(n) :: Vleft
  real, dimension(n) :: Vright

  integer, save :: lwork = -1
  integer, save :: last_n = -1

  integer :: info

  integer :: i

  ! commpute optimal workspace size only when size changes
  if(n/=last_n)then
    allocate(Work(1))
    call dgeev('N', &
               'N', &
                n, &
                A, &
                n, &
                Wreal, &
                Wimag, &
                Vleft, &
                n, &
                Vright, &
                n, &
                Work, &
                -1, &
                info)
    lwork=int(Work(1))
    last_n = n
    deallocate(Work)
  end if

  allocate(work(lwork))
  call dgeev('N', &
            'N', &
            n, &
            A, &
            n, &
            Wreal, &
            Wimag, &
            Vleft, &
            n, &
            Vright, &
            n, &
            work, &
            lwork, &
            info)

  do i=1,n
    eigenvalues(i) = COMPLEX(Wreal(i),Wimag(i))
  end do

  deallocate(work)

end subroutine

end module
