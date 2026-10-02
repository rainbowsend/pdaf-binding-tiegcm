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
! High-level driver that sets up the result writer and, each step, computes and writes the configured output quantities to NetCDF.

module result_writer_frontend

  ! intern
  use result_file_writer_module, only: nc_tiegcm, result_file_writer
  use uset_module, only: char_uset
  use configuration, only: cfg_log

  implicit none

  type(result_file_writer), target :: result_writer
  type(nc_tiegcm), target :: tiegcm_group

  logical :: t_dim_lock = .false.

  contains

  !> Registers the state and mandatory output fields as quantities to write, then creates the result writer and the main TIE-GCM grid output group.
  subroutine init_result_writer_frontend

    ! intern
    use configuration, only: cfg_output
    use state_module, only: state_vector

    use quantity_computation_module, only: qset

    implicit none

    t_dim_lock = cfg_output%lock_nc_time

    if(allocated(state_vector%tgcm_field_names))  then
      call qset%add(state_vector%tgcm_field_names, "state")
    end if

    call qset%add(cfg_output%saved_fields, "mandatory")

    call result_writer%create(trim(cfg_output%result_file_name_tag),&
                              strategy=cfg_output%output_strategy,&
                              lock_time_dim_for_locations=t_dim_lock)

    result_writer%sync_every = cfg_output%sync_every

    call tiegcm_group%init( domain_name='tiegcm_grid',&
                            save_members=cfg_output%save_members,&
                            kmax=cfg_output%max_moment,&
                            write_every_sec=cfg_output%write_every_sec,&
                            force_write_on_update=cfg_output%force_write_on_update,&
                            save_n_steps_after_update=cfg_output%save_n_steps_after_update)
    call result_writer%add(tiegcm_group)

  end subroutine

  !> Records the run's exit status and elapsed run time, then closes the result writer's output files.
  subroutine finalize_result_writer_frontend
    implicit none
    call result_writer%set_exit_status(0)
    call result_writer%write_run_time()
    call result_writer%close()
  end subroutine


  !> Declares ZG/ZGMID and all configured output quantities as NetCDF variables on every output location.
  subroutine init_result_writer

    ! extern
    use netcdf, only: NF90_COLLECTIVE, NF90_INDEPENDENT

    ! intern
    use quantity_info_module, only: LEVEL_INT,LEVEL_MID
    use quantity_computation_module, only: qset

    implicit none

    ! local
    integer :: i,j

    call tiegcm_group%add_variables("ZG")
    call tiegcm_group%add_variables("ZGMID")

    call qset%print()

    do i=1,result_writer%n_locations
      do j=1,qset%output%size()
        call result_writer%locations(i)%ptr%add_variables(qset%output%at(j))
      end do
    end do

    call result_writer%set_access(NF90_COLLECTIVE)
    call result_writer%sync()

  end subroutine

  !> Once per algorithm step, decides whether any output location is due to write, recomputes result quantities if needed, and writes the TIE-GCM grid, regular-grid, and trajectory datasets that are due.
  subroutine compute_and_write_results(algorithm_step)

    ! tie-gcm
    use init_module, only: istep

    ! intern
    use georeferenced_data_module, only: georeferenced_datasets, georeferenced_data_type, n_datasets, point, regular_grid
    use mpi_moments_module, only: allocate_and_calc_moments_mpi
    use quantity_computation_module, only:  zg_int, zg_mid, qset, compute_results
    use quantity_info_module, only: quantity_type, get_info
    use uset_module, only: str_len_char_uset
    use configuration, only: cfg_output

    implicit none

    ! arguments
    character(len=*), intent(in) :: algorithm_step ! "initial" | "forecast" | "analysis" | "advance"

    ! local

    real, allocatable, dimension(:,:,:,:) :: zg_moments
    real, allocatable, dimension(:,:,:,:) :: zgmid_moments

    integer :: i
    class(georeferenced_data_type), pointer :: georef_data

    logical :: time_to_write_any
    logical :: time_to_write



    ! Each observation can have individual write frequencies
    ! Check wether any has to be written
    time_to_write_any = tiegcm_group%time_to_write(algorithm_step)
    if((time_to_write_any.eqv..false.).and.(t_dim_lock.eqv..false.)) then
      do i=1,n_datasets
        georef_data => georeferenced_datasets(i)%ptr

        select type(georef_data)
          type is(point)
            if (georef_data%writer%time_to_write(algorithm_step)) then
              time_to_write_any = .true.
              exit
            end if
          type is(regular_grid)
            if (georef_data%writer%time_to_write(algorithm_step)) then
              time_to_write_any = .true.
              exit
            end if
        end select
      end do
    end if

    if(time_to_write_any) then
      write(*,*) 'write results ', algorithm_step, " step: ", istep
      ! in case of forecast, compute_results has been called already
      if((algorithm_step=="initial") .or.&
         (algorithm_step=="advance") .or.&
         (algorithm_step=="unconstrained_analysis") .or.&
         (algorithm_step=="analysis") ) then
        call compute_results()
      end if
      call result_writer%next_step(algorithm_step)
    end if

    if(t_dim_lock)then
      time_to_write = time_to_write_any
    else
      time_to_write = tiegcm_group%time_to_write(algorithm_step)
    end if

    if(time_to_write) then
      if(cfg_log%verbose_level>0) write(*,*) 'write TIE-GCM grid'
      call tiegcm_group%next_step(algorithm_step)

      !---- compute ZG ---------------------------

      if(tiegcm_group%save_moments)then

        call allocate_and_calc_moments_mpi(zg_int,kmax=tiegcm_group%kmax,moments=zg_moments)
        call allocate_and_calc_moments_mpi(zg_mid,kmax=tiegcm_group%kmax,moments=zgmid_moments)

        if(tiegcm_group%kmax>0) then
          call tiegcm_group%write('mean','ZG',data_tgcm_intern=zg_moments(:,:,:,1))
          call tiegcm_group%write('mean','ZGMID',data_tgcm_intern=zgmid_moments(:,:,:,1))
        end if

        if(tiegcm_group%kmax>1) then
          call tiegcm_group%write('std','ZG',data_tgcm_intern=sqrt(zg_moments(:,:,:,2)))
          call tiegcm_group%write('std','ZGMID',data_tgcm_intern=sqrt(zgmid_moments(:,:,:,2)))
        end if

        if(tiegcm_group%kmax>2) then
          call tiegcm_group%write('skew','ZG',data_tgcm_intern=zg_moments(:,:,:,3))
          call tiegcm_group%write('skew','ZGMID',data_tgcm_intern=zgmid_moments(:,:,:,3))
        end if

        if(tiegcm_group%kmax>3) then
          call tiegcm_group%write('kurt','ZG',data_tgcm_intern=zg_moments(:,:,:,4))
          call tiegcm_group%write('kurt','ZGMID',data_tgcm_intern=zgmid_moments(:,:,:,4))
        end if

        if (allocated(zg_moments)) deallocate(zg_moments)
        if (allocated(zgmid_moments)) deallocate(zgmid_moments)

      end if
      if(tiegcm_group%save_members)then
        call tiegcm_group%write('members','ZG',data_tgcm_intern=zg_int)
        call tiegcm_group%write('members','ZGMID',data_tgcm_intern=zg_mid)
      end if

      call write_nc_tiegcm(tiegcm_group, qset%output)
    end if

    do i=1,n_datasets
      georef_data => georeferenced_datasets(i)%ptr
        select type(georef_data)
          type is(point)
            if(t_dim_lock)then
              time_to_write = time_to_write_any
            else
              time_to_write = georef_data%writer%time_to_write(algorithm_step)
            end if
            if(time_to_write)then
              if(cfg_log%verbose_level>0) write(*,*) 'write ', trim(georef_data%writer%domain_name)
              call georef_data%writer%next_step(algorithm_step)
              call write_nc_trajectory(georef_data, qset%output)
            end if
          type is(regular_grid)
            if(t_dim_lock)then
              time_to_write = time_to_write_any
            else
              time_to_write = georef_data%writer%time_to_write(algorithm_step)
            end if
            if(time_to_write)then
              if(cfg_log%verbose_level>0) write(*,*) 'write ', trim(georef_data%writer%domain_name)
              call georef_data%writer%next_step(algorithm_step)
              call write_nc_reg_grid(georef_data, qset%output)
            end if
        end select

    end do

    if(time_to_write_any) then
      call result_writer%sync( (cfg_output%enforce_sync_after_update.eqv..true.) .and.&
                               (algorithm_step=="analysis") )
    end if

  end subroutine

  !> Writes the given set of quantities to the TIE-GCM grid output group, as member fields and/or ensemble moments.
  subroutine write_nc_tiegcm(spatial,set)

    ! intern
    use quantity_info_module, only: quantity_info, get_info
    use configuration, only: cfg_output
    use mpi_moments_module, only: allocate_and_calc_moments_mpi

    implicit none

    ! arguments
    type(nc_tiegcm) :: spatial
!     type(quantity_type), dimension(:), intent(in) :: quantities
    type(char_uset), intent(in) :: set

    ! local
    real, allocatable, dimension(:,:,:,:) :: moments
    type(quantity_info), pointer:: quantity

    integer :: i

    integer :: kmax

    kmax = cfg_output%max_moment

    do i=1,set%size()
      quantity => get_info(set%at(i))
      if(spatial%save_members)then
        call spatial%write('members',set%at(i),data_tgcm_intern=quantity%data)
      end if
      if(spatial%save_moments)then

        call allocate_and_calc_moments_mpi(quantity%data,kmax=kmax,moments=moments)

        if(kmax>0) call spatial%write('mean',set%at(i),data_tgcm_intern=moments(:,:,:,1))
        if(kmax>1) call spatial%write('std',set%at(i),data_tgcm_intern=sqrt(moments(:,:,:,2)))
        if(kmax>2) call spatial%write('skew',set%at(i),data_tgcm_intern=moments(:,:,:,3))
        if(kmax>3) call spatial%write('kurt',set%at(i),data_tgcm_intern=moments(:,:,:,4))

        if( allocated(moments) )deallocate(moments)

      end if
    end do

  end subroutine

  !> Writes the given set of quantities, interpolated onto a regular-grid dataset, as member fields and/or ensemble moments.
  subroutine write_nc_reg_grid(g, set)

    ! intern
    use configuration, only: cfg_output
    use georeferenced_data_module, only: regular_grid
    use mpi_moments_module, only: allocate_and_calc_moments_mpi

    implicit none

    ! arguments
    type(regular_grid) :: g
    type(char_uset), intent(in) :: set

    ! local
    real, contiguous, pointer :: dst(:,:,:)

    integer :: i

    real, allocatable, dimension(:,:,:,:) :: moments


    integer :: kmax

    kmax = cfg_output%max_moment

    do i=1,set%size()

      call g%data%get(set%at(i),dst)

      if(g%writer%save_members)then
        call g%writer%write('members',set%at(i),dst)
      end if
      if(g%writer%save_moments)then

        call allocate_and_calc_moments_mpi(dst,kmax=kmax,moments=moments)

        if(kmax>0) call g%writer%write('mean',set%at(i),moments(:,:,:,1))
        if(kmax>1) call g%writer%write('std',set%at(i),sqrt(moments(:,:,:,2)))
        if(kmax>2) call g%writer%write('skew',set%at(i),moments(:,:,:,3))
        if(kmax>3) call g%writer%write('kurt',set%at(i),moments(:,:,:,4))

        if(allocated(moments)) deallocate(moments)

      end if
    end do

  end subroutine

  !> If the trajectory dataset is valid at the current model time, writes its position and the given set of quantities as member values and/or ensemble moments.
  subroutine write_nc_trajectory(p, set)

    ! intern
    use configuration, only: cfg_output
    use georeferenced_data_module, only: point
    use quantity_info_module, only: supports_point_interpolation
    use time_module, only: get_current_modeltime
    use mpi_moments_module, only: allocate_and_calc_moments_mpi

    implicit none
    ! arguments
    type(point) :: p
    type(char_uset), intent(in) :: set

    ! local
    real, pointer :: dst
    real :: val
    real, allocatable, dimension(:) :: moments
    real :: modeltime
    logical :: is_available
    integer :: i
    real, dimension(3) :: position
    integer :: kmax

    call get_current_modeltime(modeltime)
    is_available = p%trajectory%valid_at_epoch(modeltime)

    if(is_available)then
      kmax = cfg_output%max_moment

      call p%trajectory%get_at_epoch(modeltime,val,position)
      call p%writer%write_position(position)


      do i=1,set%size()

        ! not part of a trajectory file, see nc_trajectory_add_variables
        if(supports_point_interpolation(set%at(i)).eqv..false.) cycle

        dst => null()
        call p%data%get(set%at(i),dst)

        if(p%writer%save_members)then
          call p%writer%write('members',set%at(i),dst)
        end if
        if(p%writer%save_moments)then

          call allocate_and_calc_moments_mpi(dst,kmax=kmax,moments=moments)

          if(kmax>0) call p%writer%write('mean',set%at(i),moments(1))
          if(kmax>1) call p%writer%write('std',set%at(i),sqrt(moments(2)))
          if(kmax>2) call p%writer%write('skew',set%at(i),moments(3))
          if(kmax>3) call p%writer%write('kurt',set%at(i),moments(4))

          if(allocated(moments)) deallocate(moments)

        end if
      end do

    end if
  end subroutine

end module
