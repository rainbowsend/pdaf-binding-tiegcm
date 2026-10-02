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
! Finalizes PDAF and deallocates/closes all coupling-layer modules at the end of the program.

!BOP
!
! !ROUTINE: finalize_pdaf --- Finalize PDAF
!
! !INTERFACE:
SUBROUTINE finalize_pdaf()
! !DESCRIPTION:
! This routine call MPI_finalize
!
! !USES:

  ! pdaf
  use pdaf, only: PDAF_print_info, PDAF_deallocate

  ! intern
  use configuration,&
      only: cfg_filter, deallocate_configuration
  use mod_parallel_pdaf,&
      only: rank_world
  use model_parameter_IO_module,&
      only: model_parameter_writer
  use model_parameter_handling_module,&
      only: deallocate_model_parameter_handling_module
  use observations_module,&
      only: deallocate_observations_module
  use quantity_computation_module,&
      only: deallocate_quantity_computation_module
  use quantity_info_module,&
      only: finalize_quantity_info_module
  use result_writer_frontend,&
      only: finalize_result_writer_frontend
  use state_module,&
      only: deallocate_state_module
  use structured_gird_subdomain_module,&
      only: deallocate_structured_gird_subdomain_module
  use time_module, only: finalize_time_module

  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: tgcm.F program
!EOP

! in open loop no PDAF intern variables are allocated
  if(cfg_filter%open_loop .eqv. .false.) then
  ! *** Show allocated memory for PDAF ***
    IF (rank_world==0) CALL PDAF_print_info(2)

  ! *** Print PDAF timings onto screen ***
    IF (rank_world==0) CALL PDAF_print_info(3)

  ! *** Deallocate PDAF arrays
    CALL PDAF_deallocate()
  end if

  call deallocate_state_module

  call deallocate_observations_module

  call deallocate_configuration

  call finalize_result_writer_frontend

  call deallocate_structured_gird_subdomain_module

  call deallocate_quantity_computation_module

  call finalize_quantity_info_module

  call deallocate_model_parameter_handling_module

  call model_parameter_writer%close()

  call finalize_time_module

END SUBROUTINE finalize_pdaf
