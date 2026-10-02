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
! Instantiates and manages the configured set of observation-type objects used for assimilation.

module observations_module

  ! intern
  use tgcm_pdaf_omi_obs_type_module, only: tgcm_pdaf_observation_interface
  use obs_den_grid_pdafomi, only: obs_den_grid
  use obs_tum_ne_pdafomi, only: obs_tum_ne
  use obs_tum_vtec_pdafomi, only: obs_tum_vtec
  use obs_satellite_pdafomi, only: obs_satellite
  use configuration, only: max_cal_den_tme_file, max_satellite_file, cfg_filter

  implicit none

  type :: observation
    class(tgcm_pdaf_observation_interface), pointer :: ptr
  end type

  ! contains all initalized observations
  type(observation), dimension(:), allocatable :: observations
  ! subset of observations array, that are assimilated
  type(observation), dimension(:), allocatable :: observations_assim

  type(obs_den_grid), dimension(max_cal_den_tme_file), target :: msis_cal

  type(obs_tum_ne), target :: tum_ne

  type(obs_tum_vtec), target :: tum_vtec

  type(obs_satellite), dimension(max_satellite_file), target :: satellite

  integer :: n_assim_obs

  contains

  !> Allocates and initializes all configured observation-type instances
  !! (msis_cal, tum_ne, satellite) and builds the observations/
  !! observations_assim arrays.
  subroutine init_observations_module

    ! intern
    use configuration, only: cfg_observation

    implicit none

    ! local
    integer :: n
    integer :: i, obs_type_idx, msis_cal_idx, satellite_idx

    logical :: assim_run

    logical, dimension(:), allocatable :: mask

    assim_run = cfg_filter%open_loop.eqv..false.

    if(cfg_filter%open_loop) then
      n = count(cfg_observation%cal_den%always_save) + &
          count(cfg_observation%satellite%always_save) + &
          count((/cfg_observation%tum_ne%always_save/)) + &
          count((/cfg_observation%tum_vtec%always_save/))
    else
      n = count((cfg_observation%cal_den%apply .or. cfg_observation%cal_den%always_save)) + &
          count((cfg_observation%satellite%apply .or. cfg_observation%satellite%always_save)) + &
          count((/(cfg_observation%tum_ne%apply .or. cfg_observation%tum_ne%always_save)/)) + &
          count((/(cfg_observation%tum_vtec%apply .or. cfg_observation%tum_vtec%always_save)/))
    end if

    if(cfg_filter%open_loop) then
      n_assim_obs = 0
    else
      n_assim_obs = count(cfg_observation%cal_den%apply) + &
                    count(cfg_observation%satellite%apply) + &
                    count((/cfg_observation%tum_ne%apply/)) + &
                    count((/cfg_observation%tum_vtec%apply/))
    end if

    write(*,'(i2,a,i2,a)') n, ' observation types are initalized of which ', n_assim_obs , ' are assimilated'

    allocate(observations(n))
    allocate(mask(n))
    mask = .false.

    obs_type_idx=0
    msis_cal_idx=0
    do i =1,size(cfg_observation%cal_den)
      if((cfg_observation%cal_den(i)%apply .and. assim_run) .or. &
         (cfg_observation%cal_den(i)%always_save)) then

        obs_type_idx=obs_type_idx+1
        msis_cal_idx=msis_cal_idx+1
        call msis_cal(msis_cal_idx)%init(cfg_observation%cal_den(i))
        observations(obs_type_idx)%ptr => msis_cal(msis_cal_idx)
        mask(obs_type_idx) = cfg_observation%cal_den(i)%apply

      end if
    end do

    if((cfg_observation%tum_ne%apply .and. assim_run) .or. &
       (cfg_observation%tum_ne%always_save)) then

      obs_type_idx=obs_type_idx+1
      call tum_ne%init(cfg_observation%tum_ne)
      observations(obs_type_idx)%ptr => tum_ne
      mask(obs_type_idx) = cfg_observation%tum_ne%apply

    end if

    if((cfg_observation%tum_vtec%apply .and. assim_run) .or. &
       (cfg_observation%tum_vtec%always_save)) then

      obs_type_idx=obs_type_idx+1
      call tum_vtec%init(cfg_observation%tum_vtec)
      observations(obs_type_idx)%ptr => tum_vtec
      mask(obs_type_idx) = cfg_observation%tum_vtec%apply

    end if

    satellite_idx = 0
    do i =1,size(cfg_observation%satellite)
      if((cfg_observation%satellite(i)%apply .and. assim_run) .or. &
         (cfg_observation%satellite(i)%always_save)) then

        obs_type_idx=obs_type_idx+1
        satellite_idx=satellite_idx+1
        call satellite(satellite_idx)%init(cfg_observation%satellite(i))
        observations(obs_type_idx)%ptr => satellite(satellite_idx)
        mask(obs_type_idx) = cfg_observation%satellite(i)%apply

      end if
    end do

    ! every initialized observation type gets its own slot in the
    ! thread-private local_obs_pool, see pdaf_omi_obs_type_module
    do i = 1, n
      call observations(i)%ptr%set_obs_l_id(i)
    end do

    if(assim_run) then
      allocate(observations_assim(n_assim_obs),source=pack(observations,mask))
    end if

    deallocate(mask)

  end subroutine

  !> Deallocates all initialized observation-type instances and the
  !! observations/observations_assim arrays.
  subroutine deallocate_observations_module

    implicit none

    integer :: i

    if(allocated(observations)) then
      do i=1, size(observations)
        call observations(i)%ptr%deallocate
      end do
      deallocate(observations)
      if(allocated(observations_assim)) deallocate(observations_assim)
    end if

  end subroutine

end module observations_module
