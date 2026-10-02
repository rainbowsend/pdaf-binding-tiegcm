! Based on PDAF template/tutorial code.
! Copyright (c) 2004-2026 Lars Nerger, Alfred Wegener Institute,
! Helmholtz Center for Polar and Marine Research, Bremerhaven, Germany.
!
! Modified for the TIE-GCM/PDAF coupling.
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
! PDAF callback: sets the number of local analysis domains (subdomains) on this process.

!BOP
!
! !ROUTINE: init_n_domains_pdaf --- Set number of local analysis domains
!
! !INTERFACE:
SUBROUTINE init_n_domains_pdaf(step, n_domains_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called in PDAF\_X\_update 
! at the beginning of the analysis step before 
! the loop through all local analysis domains. 
! It has to set the number of local analysis 
! domains for the PE-local domain.
!
! The routine is called by all filter processes.
!
! !USES:

  ! intern
  use structured_gird_subdomain_module, only: nsub

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: step        ! Current time step
  INTEGER, INTENT(out) :: n_domains_p ! PE-local number of analysis domains

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_n_domains)
! Called by: PDAF_lestkf_update  (as U_init_n_domains)
! Called by: PDAF_letkf_update   (as U_init_n_domains)
! Called by: PDAF_lnetf_update   (as U_init_n_domains)
!EOP


! ************************************
! *** Initialize number of domains ***
! ************************************

  ! Only state is localized. Global parameters are updated later in prepoststep

  n_domains_p = nsub

END SUBROUTINE init_n_domains_pdaf
