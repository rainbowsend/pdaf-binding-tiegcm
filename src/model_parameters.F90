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
! Initializes model parameters from ensemble files at startup and applies perturbations to TIE-GCM input constants.

module model_parameter_module

use model_parameter_handling_module, only: model_parameter

implicit none

contains

!> Reads the ensemble of model parameters from their configured NetCDF files, sets the TIE-GCM start time, and applies each parameter's initial handling (perturb/calibrate/overwrite).
subroutine init_parameter_module

  ! tie-gcm
  use input_module,only: pristart
  use hist_module,only: modeltime

  ! intern
  use model_parameter_handling_module, only:&
    init_model_parameter_handling_module, &
    model_parameters,&
    PARAM_HANDLING_NONE,&
    PARAM_HANDLING_PERTURB,&
    PARAM_HANDLING_CALIBRATE,&
    PARAM_HANDLING_OVERWRITE,&
    write_parmeter_handling,&
    n_parameters
  use model_parameter_IO_module, only: model_parameter_reader_type
  use configuration, only: cfg_parameters

  implicit none

  ! local
  type(model_parameter_reader_type) :: reader_default
  type(model_parameter_reader_type) :: reader

  integer :: i
  integer, dimension(:), allocatable :: skip_list

  call init_model_parameter_handling_module

  call get_unique_sorted_skip_list(skip_list)
  if(size(skip_list)>0)then
    write(*,*) "ATTENTION the perturbations with the following index are skipped:  ", skip_list
  end if

  write(*,*) 'reading ensemble of model parameters (dynamics)'

  call reader_default%open(cfg_parameters%ensemble_file, skip_list)
  do i=1,n_parameters
    if (model_parameters(i)%ensemble_file == "default") then
        write(*,*) 'read ', model_parameters(i)%name, ' from ', cfg_parameters%ensemble_file
        call reader_default%read(model_parameters(i))
    else
        write(*,*) 'read ', model_parameters(i)%name, ' from ', model_parameters(i)%ensemble_file
        call reader%open(model_parameters(i)%ensemble_file, skip_list)
        call reader%read(model_parameters(i))
        call reader%close
    end if
  end do

  call reader_default%close

  ! modeltime is set in advacne.F for the first time. However we need it already to be set here
  modeltime = 0
  modeltime(1:4) = pristart(:,1)

  do i=1, n_parameters
    call model_parameters(i)%update()
  end do

  call write_parmeter_handling()

  if(allocated(skip_list)) deallocate(skip_list)


end subroutine

!> Builds the list of indices of pertubrations in the NetCDF perturbation file that are skipped. Only works if the perturbation file has more members than the model ensemble.
subroutine get_unique_sorted_skip_list(list)

  use m_unirnk, only: unirnk

  use configuration, only: cfg_parameters, skip_list_len

  implicit none

  ! args
  integer, dimension(:), allocatable, intent(inout) :: list

  ! local
  integer, dimension(:), allocatable :: valid, unique_id
  logical, dimension(skip_list_len) :: mask
  integer :: n

  mask = .false.
  where(cfg_parameters%skip_list>0)
    mask = .true.
  end where
  allocate(valid(count(mask)))
  allocate(unique_id(count(mask)))
  valid = pack(cfg_parameters%skip_list,mask)

  call unirnk(valid, unique_id, n)

  if(allocated(list)) deallocate(list)
  allocate(list(n))
  list=valid(unique_id(1:n))

  deallocate(valid,unique_id)

end subroutine

!> Resets and reapplies the configured perturbation to the TIE-GCM collision/Joule-heating factor.
subroutine apply_perturbation_to_input_constants

  ! intern
  use model_parameter_handling_module, only: &
    param_colfac,param_joulefac, &
    PARAM_HANDLING_PERTURB

  ! tie-gcm
  use input_module, only: colfac, colfac_0, joulefac, joulefac_0

  implicit none

  if(param_colfac%handling==PARAM_HANDLING_PERTURB) colfac = colfac_0
  call param_colfac%handle(colfac)

  if(param_joulefac%handling==PARAM_HANDLING_PERTURB) joulefac = joulefac_0
  call param_joulefac%handle(joulefac)

end subroutine apply_perturbation_to_input_constants


end module model_parameter_module
