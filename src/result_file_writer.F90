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
! Core parallel NetCDF output engine: defines result-file/spatial-domain types and writes TIE-GCM grid, regular-grid, and trajectory data (incl. ensemble moments) to disk.

! strategy (1) collect everything on world root and write (there is no parallel file system, e.g. NFS)
!          (2) use netcdf parallel IO feature (there is a parallel file system e.g. GPFS). This is currently not faster than strategy 1. do not use
!
! each quantity/variable (e.g., neutral density) is associated with a
! spatial domain (e.g., TIE-GCM grid or gridded MSIS observations)
! and with one of the following types: 'mean', 'std', 'skew', 'kurt' or 'members'.
! We use groups in the netcdf file to represent this. An example of the
! hierachy is given below
!                                           |
!                       -----------------------------------------
!                      /                                         \
!                     /                                           \
!           ----TIE-GCM-grid -----                       ----  obs-grid  -----
!          /          |           \                     /          |           \
!         /           |            \                   /           |            \
!        /            |             \                 /            |             \
!      mean          std          members           mean          std          members
!      /|\           /|\            /|\             /|\           /|\            /|\
!     / | \         / | \          / | \           / | \         / | \          / | \
!    /  |  \       /  |  \        /  |  \         /  |  \       /  |  \        /  |  \
!  den  ne  tn   den  ne  tn    den  ne  tn     den  ne  tn   den  ne  tn    den  ne  tn
!
! The leave gropus (mean, std, members) can be located on different NetCDF files.
! For example, one file could include all members and a second file all
! moments (mean and std):
!
!              /                   |                   |                  \
!             /                    |                   |                   \
!       TIE-GCM-grid          obs-grid-1          obs-grid-2          obs-trajectory-n
!            /|\                  /|\                 /|\                  /|\
!           / | \                / | \               / | \                / | \
!          /  |  \              /  |  \             /  |  \              /  |  \
!      mean  std  members   mean  std  members  mean  std  members   mean  std  members
!        |    |     |         |    |     |        |    |     |         |    |     |
!        v    v     V         v    v     V        v    v     V         v    v     V
!     file1 file1 file2    file1 file1 file2   file1 file1 file2    file3 file3 file3
!
!
! Short introduction to this module
!
! 1) You need an instance of result_file_writer that can be assigned to a pointer:
!       type(result_file_writer), target :: rfw
!
! 2) You need one or more instances of a spatial domain
!       type(nc_tiegcm) :: tiegcm_group
!
! 3) Add one or more files to the result_file_writer
!       call rfw%add_file('file1.nc')
!       call rfw%add_file('file2.nc')
!
! 4) link the spatial domain to the writer
!       call rfw%link(tiegcm_group)
!
! 5) assign to each spatial domain and each type(mean, std, members)
!    one of the files created in step3
!       call tiegcm_group%assign(rfw%files(1),"mean")
!       call tiegcm_group%assign(rfw%files(1),"std")
!       call tiegcm_group%assign(rfw%files(2),"members")
!
! 6) at each iteration increment the temporal dimension
!       call rfw%next_step('advance')
!
! 7) write data calling the write subroutine
!
!       call tiegcm_group%write('mean','ZG',data_tgcm_intern=testdata)
!
! 8) at the end of the program close all files
!       call rfw%close()
!
module result_file_writer_module

  ! extern
  use netcdf

  ! tie-gcm
  use nchist_module, only: handle_ncerr

  ! intern
  use trajectory_data_module, only: trajectory_data
  use grid_observation_module, only: reg_grid_dataset_root

  use field_bundle_module, only: bundle

  implicit none

  integer, parameter :: SINGLE_FILE_SINGLE_WRITER = 1
  integer, parameter :: MULTIPLE_FILE_MULTIPLE_WRITER = 2

  character(len=4), dimension(4), parameter :: moment_nc_name = (/"mean","std ","skew","kurt"/)

  type, abstract :: spatial_domain
    ! member variables
    character(len=64) :: domain_name
    logical :: save_members = .false.
    logical :: save_moments = .false.
    integer :: kmax = 2
    integer, private :: root_id_members = -1
    integer, private, dimension(4) :: root_id_moments = -1
    integer, private :: spatial_id_members = -1
    integer, private, dimension(4) :: spatial_id_moments = -1
    integer, private :: grp_id_members = -1
    integer, private, dimension(4) :: grp_id_moments = -1
    type(result_file_writer), pointer, private :: result_writer
    integer :: counter = 0
    integer, private :: write_every = -1
    logical, private :: force_write_on_update = .true.
    integer, private :: save_n_steps_after_update = 0
    integer, private :: istep_last_update = -1
    integer :: istep_last_write = -1
    contains
    procedure, pass(this) :: assign => spatial_domain_assign
    procedure, pass(this) :: init => spatial_domain_init
    procedure, pass(this) :: get_grp_id => spatial_domain_get_grp_id
    procedure, pass(this) :: get_unique_spatial_ids => spatial_domain_get_unique_spatial_ids
    procedure, pass(this) :: next_step => spatial_domain_next_step
    procedure, pass(this) :: get_counter => spatial_domain_get_counter
    procedure, pass(this) :: time_to_write => spatial_domain_time_to_write
    procedure(add_dim_interface), deferred :: add_dimensions
    procedure(add_meta_data_interface), deferred :: add_meta_data
    procedure(add_var_interface), deferred :: add_variables
  end type

  abstract interface
    !> Deferred interface: adds a spatial domain's NetCDF spatial dimensions to a group.
    subroutine add_dim_interface(this, ncid)
          import spatial_domain
          class(spatial_domain), intent(inout) :: this
          integer, intent(in) :: ncid
    end subroutine
    !> Deferred interface: adds a spatial domain's NetCDF metadata attributes to a group.
    subroutine add_meta_data_interface(this, ncid)
          import spatial_domain
          class(spatial_domain), intent(inout) :: this
          integer, intent(in) :: ncid
    end subroutine
    !> Deferred interface: defines one NetCDF variable for a spatial domain.
    subroutine add_var_interface(this,varname)
          import spatial_domain
          class(spatial_domain) :: this
          character(len=*), intent(in) :: varname
    end subroutine
  end interface

  type, extends(spatial_domain) :: nc_tiegcm
    contains
    procedure, pass(this) :: add_dimensions => nc_tiegcm_add_dimensions
    procedure, pass(this) :: add_variables => nc_tiegcm_add_variables
    procedure, pass(this) :: add_meta_data => nc_tiegcm_add_meta_data
    procedure, pass(this) :: write => nc_tiegcm_write
  end type nc_tiegcm

  type, extends(spatial_domain) :: nc_reg_grid
    type(reg_grid_dataset_root), pointer :: dataset
    contains
    procedure, pass(this) :: add_dimensions => nc_reg_grid_add_dimensions
    procedure, pass(this) :: add_variables => nc_reg_grid_add_variables
    procedure, pass(this) :: add_meta_data => nc_reg_grid_add_meta_data
    procedure, pass(this) :: write => nc_reg_grid_write
    procedure, pass(this) :: link_dataset => nc_reg_grid_link_dataset
  end type nc_reg_grid

  type, extends(spatial_domain) :: nc_trajectory
    type(trajectory_data), pointer :: dataset
    contains
    procedure, pass(this) :: add_dimensions => nc_trajectory_add_dimensions
    procedure, pass(this) :: add_variables => nc_trajectory_add_variables
    procedure, pass(this) :: add_meta_data => nc_trajectory_add_meta_data
    procedure, pass(this) :: write => nc_trajectory_write
    procedure, pass(this) :: write_position => nc_trajectory_write_position
    procedure, pass(this) :: link_dataset => nc_trajectory_link_dataset
  end type nc_trajectory

  type :: spatial
    class(spatial_domain), pointer :: ptr
  end type

  type :: result_file
    integer :: ncid = -1
    character(len=128) :: name
    logical :: parallel_access
    logical :: lock_time_dim_for_locations
    contains
    procedure, pass(this) :: create => result_file_create
    procedure, pass(this) :: close => result_file_close
    procedure, pass(this) :: next_step => result_file_next_step
    procedure, pass(this) :: sync => result_file_sync
    procedure, pass(this) :: set_exit_status => result_file_set_exit_status
    procedure, pass(this) :: write_run_time => result_file_write_run_time
  end type

  integer, parameter :: max_files = 3
  integer, parameter :: max_locations = 10

  type :: result_file_writer
    type(result_file), dimension(max_files) :: files
    integer :: n_files = 0
    !
    type(spatial), dimension(max_locations) :: locations
    integer :: n_locations = 0
    !
    integer :: counter = 0
    !
    integer :: strategy
    !
    integer :: sync_every = 1
    !
    logical :: lock_time_dim_for_locations
    contains

    procedure, pass(this) :: create => result_file_writer_create
    procedure, pass(this) :: add_file => result_file_writer_add_file
    procedure, pass(this) :: add => result_file_writer_add
    procedure, pass(this) :: close => result_file_writer_close
    procedure, pass(this) :: next_step => result_file_writer_next_step
    procedure, pass(this) :: sync => result_file_writer_sync
    procedure, pass(this) :: link => result_file_writer_link
    procedure, pass(this) :: set_exit_status => result_file_writer_set_exit_status
    procedure, pass(this) :: write_run_time => result_file_writer_write_run_time
    procedure, pass(this) :: set_access => result_file_writer_set_access
  end type

  integer, protected :: nc_var_xtype = NF90_FLOAT ! NF90_FLOAT or NF90_DOUBLE
  logical, parameter :: nc_var_shuffle = .true.
  logical, parameter :: nc_var_fletcher32 = .false.
  integer, parameter :: nc_var_deflate_level = 4 ! ATTENTION values larger than one causes HDF error when performing independent parallel I/O using unlimited dimensions

contains

  !> Standalone smoke test exercising the result_file_writer/nc_tiegcm API end to end.
  subroutine test_result_file_writer_module

    ! tie-gcm
    use fields_module,only: levd0,levd1,lond0,lond1,latd0,latd1

    ! intern
    use mod_parallel_pdaf, only: rank_world

    implicit none

    ! local
    real, dimension(levd0:levd1,lond0:lond1,latd0:latd1):: testdata

    type(nc_tiegcm) :: tiegcm_group

    type(result_file_writer), target :: rfw

    call rfw%add_file('file1.nc')
    call rfw%add_file('file2.nc')

    call tiegcm_group%init('tiegcm_grid',.true.,4,1,.true.,0)

    call rfw%link(tiegcm_group)

    if(rank_world==0) then
      call tiegcm_group%assign(rfw%files(1),"mean")
      call tiegcm_group%assign(rfw%files(1),"std")
      call tiegcm_group%assign(rfw%files(1),"skew")
      call tiegcm_group%assign(rfw%files(1),"kurt")
      call tiegcm_group%assign(rfw%files(2),"members")
    end if

    call tiegcm_group%add_variables("ZG")

    call rfw%next_step('advance')

    testdata=rank_world
    call tiegcm_group%write('mean','ZG',data_tgcm_intern=testdata)
    call tiegcm_group%write('members','ZG',data_tgcm_intern=testdata)

    call rfw%close()

  end subroutine

  !> Creates the writer's output file(s) according to the chosen write strategy (single vs. multiple files/writers).
  subroutine result_file_writer_create(this, basename, strategy, lock_time_dim_for_locations)

    ! extern
    use mpi_f08

    ! intern
    use mod_parallel_pdaf, only: COMM_filter, rank_world, filterpe

    implicit none

! required for nf_set_log_level
!#include <netcdf.inc>

    ! arguments
    class(result_file_writer) :: this
    character(len=*), intent(in) :: basename
    integer, intent(in) :: strategy
    logical, intent(in) :: lock_time_dim_for_locations

    ! local
    character(len=128) :: error_message

!     integer :: n
!     integer :: istat

    ! netcdf library must be compiled with '--enable-logging' for that feature
!     istat = nf_set_log_level(0)

    this%strategy = strategy

    this%lock_time_dim_for_locations = lock_time_dim_for_locations

    select case(this%strategy)
      case(SINGLE_FILE_SINGLE_WRITER)
        write(*,*) 'Using single file single writer strategy'
        ! netcdf file is opened only on root of MPI_COMM_WORLD
        if(rank_world==0) then
          call this%add_file(trim(basename)//'.nc')
        end if
      case(MULTIPLE_FILE_MULTIPLE_WRITER)
          write(*,*) 'Using two files multiple writer strategy'

          ! ATTENTION: compute_temporal_dimension_size() must be revised before
          ! the fixed time dimension (dim_t=n) below is reactivated. It sizes the
          ! time dimension as
          !     n = n_assim_step*2 + n_advance - n_skip
          ! which only accounts for cfg_output%write_every_sec and for the two
          ! update steps (forecast + analysis) of each assimilation cycle. It does
          ! NOT account for the debugging options in cfg_output:
          !   - save_n_steps_after_update adds up to
          !     n_assim_step*save_n_steps_after_update extra advance writes
          !     (see spatial_domain_time_to_write)
          !   - save_unconstrained_analysis adds a THIRD update step per cycle
          !     (see distribute_state_pdaf.F90, compute_and_write_results
          !     ("unconstrained_analysis")), which the factor 2 above does not cover
          ! With either option enabled the fixed dimension would be too small and
          ! the writes would run past the end of the time dimension. The unlimited
          ! time dimension currently in use is not affected.
!         n=compute_temporal_dimension_size()

        ! netcdf file is only opened on filter processes
        call this%add_file(trim(basename)//'.nc',&
                          !dim_t=n,&
                          comm=COMM_filter,&
                          rank_mask=filterpe)

        ! netcdf file is opened on all processes
        call this%add_file(trim(basename)//'_members.nc',&
                          ! dim_t=n,&
                           comm=MPI_COMM_WORLD)

        ! netcdf file is opened only on root all processes
        call this%add_file(trim(basename)//'_tra.nc',&
                          ! dim_t=n,&
                          rank_mask=(rank_world==0))


      case default
        write(error_message,*) "creating result file: unknown/invalid write strategy ", this%strategy
        call shutdown(trim(error_message))
    end select

  end subroutine

  !> Links a spatial domain to the writer and assigns it to the mean/std/skew/kurt/members file(s) per strategy.
  subroutine result_file_writer_add(this, spatial)

    ! extern
    use mpi_f08

    implicit none

    ! arguments
    class(result_file_writer) :: this
    class(spatial_domain), intent(inout) :: spatial

    ! local
    character(len=128) :: error_message

    call this%link(spatial)

    select case(this%strategy)
      case(SINGLE_FILE_SINGLE_WRITER)
          if(spatial%save_moments) call spatial%assign(this%files(1), "mean")
          if(spatial%save_moments) call spatial%assign(this%files(1), "std")
          if(spatial%save_moments) call spatial%assign(this%files(1), "skew")
          if(spatial%save_moments) call spatial%assign(this%files(1), "kurt")
          if(spatial%save_members) call spatial%assign(this%files(1), "members")

      case(MULTIPLE_FILE_MULTIPLE_WRITER)

          select type(spatial)
          type is(nc_reg_grid)
            if(spatial%save_moments) call spatial%assign(this%files(1), "mean")
            if(spatial%save_moments) call spatial%assign(this%files(1), "std")
            if(spatial%save_moments) call spatial%assign(this%files(1), "skew")
            if(spatial%save_moments) call spatial%assign(this%files(1), "kurt")
            if(spatial%save_members) call spatial%assign(this%files(2), "members")
          type is(nc_tiegcm)
            if(spatial%save_moments) call spatial%assign(this%files(1), "mean")
            if(spatial%save_moments) call spatial%assign(this%files(1), "std")
            if(spatial%save_moments) call spatial%assign(this%files(1), "skew")
            if(spatial%save_moments) call spatial%assign(this%files(1), "kurt")
            if(spatial%save_members) call spatial%assign(this%files(2), "members")
          type is(nc_trajectory)
            if(spatial%save_moments) call spatial%assign(this%files(3), "mean")
            if(spatial%save_moments) call spatial%assign(this%files(3), "std")
            if(spatial%save_moments) call spatial%assign(this%files(3), "skew")
            if(spatial%save_moments) call spatial%assign(this%files(3), "kurt")
            if(spatial%save_members) call spatial%assign(this%files(3), "members")
          end select

      case default
        write(error_message,*) "adding file to result file: unknown/invalid write strategy ", this%strategy
        call shutdown(trim(error_message))
    end select

  end subroutine

  !> Creates and appends one more result_file, optionally scoped to a given MPI communicator/rank subset.
  subroutine result_file_writer_add_file(this, filename, comm, dim_t, rank_mask)

    ! extern
    use mpi_f08

    implicit none

    ! arguments
    class(result_file_writer) :: this
    character(len=*), intent(in) :: filename
    type(MPI_Comm), intent(in), optional :: comm
    integer, intent(in), optional :: dim_t
    logical, intent(in), optional :: rank_mask ! only create files on rank where this argument is true

    ! local
    integer :: dim_t_

    if(present(dim_t)) then
      dim_t_ = dim_t
    else
      dim_t_ = -1
    end if

    this%n_files = this%n_files+1

    if(present(rank_mask)) then
      if(rank_mask .eqv. .false.) return
    end if

    if(present(comm)) then
      call this%files(this%n_files)%create(trim(filename), dim_t_, this%lock_time_dim_for_locations, comm)
    else
      call this%files(this%n_files)%create(trim(filename), dim_t_, this%lock_time_dim_for_locations)
    end if

  end subroutine

  !> Registers a spatial domain as one of the writer's output locations.
  subroutine result_file_writer_link(this, spatial_group)

    implicit none

    ! arguments
    class(result_file_writer), target :: this
    class(spatial_domain), target :: spatial_group

    ! local
    character(len=128) :: error_message

    spatial_group%result_writer => this

    this%n_locations = this%n_locations + 1

    if(this%n_locations <= max_locations) then
      this%locations(this%n_locations)%ptr => spatial_group
    else
      write(error_message,'(a,i3,a)') 'result_file_writer_link: exceeding maximal number of locations (',&
        max_locations, '). Increase it accordingly.'
      call shutdown(trim(error_message))
    end if

  end subroutine

  !> Closes all result files owned by the writer.
  subroutine result_file_writer_close(this)

    implicit none

    ! arguments
    class(result_file_writer) :: this

    ! local
    integer :: i

    do i = 1, this%n_files
      call this%files(i)%close()
    end do

  end subroutine

  !> Writes the run's exit-status flag to all result files.
  subroutine result_file_writer_set_exit_status(this,status)

    implicit none
    ! arguments
    class(result_file_writer) :: this
    integer, intent(in) :: status

    ! local
    integer :: i

    do i = 1, this%n_files
      call this%files(i)%set_exit_status(status)
    end do

  end subroutine

  !> Writes the accumulated total run time to all result files.
  subroutine result_file_writer_write_run_time(this)

    implicit none
    ! arguments
    class(result_file_writer) :: this

    ! local
    integer :: i

    do i = 1, this%n_files
      call this%files(i)%write_run_time()
    end do

  end subroutine

  !> Advances the shared temporal-dimension counter and records it in all result files.
  subroutine result_file_writer_next_step(this, algorithm_step)

    implicit none

    ! arguments
    class(result_file_writer) :: this
    character(len=*), intent(in) :: algorithm_step

    ! local
    integer :: i

    if(this%lock_time_dim_for_locations) then

      this%counter = this%counter+1

      do i = 1, this%n_files
        call this%files(i)%next_step(algorithm_step, this%counter)
      end do
    end if

  end subroutine

  !> Flushes (syncs) all result files to disk, throttled by sync_every.
  subroutine result_file_writer_sync(this,enforce)

    implicit none

    ! arguments
    class(result_file_writer) :: this
    logical, optional, intent(in) :: enforce

    ! local
    integer :: i
    integer :: counter
    integer, save :: counter_last = -1

    logical :: enforce_

    if(present(enforce)) then
      enforce_ = enforce
    else
      enforce_ = .false.
    end if

    if(this%lock_time_dim_for_locations) then
      counter = this%counter
    else
      ! ATTENTION currently tiegcm_grid counter is used
      counter = this%locations(1)%ptr%counter
    end if

    if((MODULO(counter,this%sync_every)==0).and.(counter_last/=counter)) then
      write(*,*) 'syncing nc files ', counter
      do i = 1, this%n_files
        call this%files(i)%sync
      end do
      counter_last = counter
    else if(enforce_) then
      write(*,*) 'syncing nc files was enforced'
      do i = 1, this%n_files
        call this%files(i)%sync
      end do
    end if
  end subroutine

  !> Switches NetCDF variables between collective and independent parallel-I/O access.
  subroutine result_file_writer_set_access(this,var_access)

    implicit none

    ! arguments
    class(result_file_writer) :: this
    integer, intent(in) :: var_access

    ! local
    integer :: i
    integer :: ncid

    if(this%strategy==MULTIPLE_FILE_MULTIPLE_WRITER) then

      ! temporal coordinate dimension vars are all saved in root group
      !
      ! They must be collective in any case (collective or independent). Else we get:
      ! 'Attempt to extend dataset during NC_INDEPENDENT I/O operation. Use nc_va'
      do i = 1, this%n_files
        if(this%files(i)%parallel_access) then
          call set_collective_access_of_all_vars(this%files(i)%ncid)
        end if
      end do

      if(var_access==NF90_COLLECTIVE)then
          do i=1,this%n_locations
            ncid = this%locations(i)%ptr%grp_id_members
            call set_collective_access_of_all_vars(ncid)
          end do

          do i=1,this%n_locations
            ncid = this%locations(i)%ptr%grp_id_moments(1)
            call set_collective_access_of_all_vars(ncid)
          end do

          do i=1,this%n_locations
            ncid = this%locations(i)%ptr%grp_id_moments(2)
            call set_collective_access_of_all_vars(ncid)
          end do

          do i=1,this%n_locations
            ncid = this%locations(i)%ptr%grp_id_moments(3)
            call set_collective_access_of_all_vars(ncid)
          end do

          do i=1,this%n_locations
            ncid = this%locations(i)%ptr%grp_id_moments(4)
            call set_collective_access_of_all_vars(ncid)
          end do
      end if
    end if

  end subroutine

  !> Sets NF90_COLLECTIVE parallel access mode on every variable in a NetCDF group.
  subroutine set_collective_access_of_all_vars(ncid)
    implicit none

    ! arguments
    integer, intent(in) :: ncid

    ! local
    integer :: i
    integer :: nvars
    integer, dimension(128) :: var_buff
    integer :: istat

    if(ncid >= 0) then
      istat = nf90_inq_varids(ncid,nvars,var_buff)
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'set_collective_access_of_all_vars inq')
      do i=1,nvars
        istat = nf90_var_par_access(ncid,var_buff(i),NF90_COLLECTIVE)
!         if (istat /= NF90_NOERR) call handle_ncerr(istat,'set_collective_access_of_all_vars')
      end do
    end if

  end subroutine

  !> Selects whether output variables are written as NF90_FLOAT or NF90_DOUBLE.
  subroutine result_file_writer_module_set_precision(use_double)

    implicit none

    ! arguments
    logical, intent(in) :: use_double

    if(use_double) then
      nc_var_xtype = NF90_DOUBLE
    else
      nc_var_xtype = NF90_FLOAT
    end if
  end subroutine

  !> Writes the algorithm-step flag (initial/forecast/analysis/...) and advances model time for one NetCDF file/group.
  subroutine next_step(ncid, algorithm_step, counter)

    ! intern
    use netcdf_functionality, only: add_model_time

    implicit none

    ! arguments
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: algorithm_step
    integer, intent(in) :: counter

    ! local
    integer :: var_id
    integer :: istat

    istat = nf90_inq_varid(ncid, "step", var_id)
    select case (algorithm_step)
      case ("initial")
        istat = nf90_put_var(ncid=ncid, varid=var_id, values=-1, start=(/counter/))
      case ("advance")
        istat = nf90_put_var(ncid=ncid, varid=var_id, values=0, start=(/counter/))
      case ("forecast")
        istat = nf90_put_var(ncid=ncid, varid=var_id, values=1, start=(/counter/))
      case ("analysis")
        istat = nf90_put_var(ncid=ncid, varid=var_id, values=2, start=(/counter/))
      case ("unconstrained_analysis")
        istat = nf90_put_var(ncid=ncid, varid=var_id, values=3, start=(/counter/))
      case default
        call shutdown("invalid input in save_fields. Use 'initial', 'forecast', 'analysis', 'unconstrained_analysis', or 'advance'")
    end select
    call add_model_time(ncid, counter, .false.)

  end subroutine

  !> Creates one NetCDF output file with its temporal/ensemble dimensions and global run metadata.
  subroutine result_file_create(this, filename, dim_t, lock_time_dim_for_locations, comm )

    ! extern
    use mpi_f08

    ! intern
    use mod_parallel_pdaf, only: n_modeltasks
    use netcdf_functionality, only: add_global_meta_data, init_temporal_dims

    implicit none

    ! arguments
    class(result_file) :: this
    character(len=*), intent(in) :: filename
    integer, intent(in) :: dim_t
    logical, intent(in) :: lock_time_dim_for_locations
    type(MPI_Comm), intent(in), optional :: comm

    ! local
    integer :: istat

    integer :: dim_id_ulim
    integer :: var_id_info
    integer :: var_id_status
    integer :: var_id_runtime
    integer :: dim_id_ens

    character(len=16384) :: file_content_buff
    character(len=1024) :: arg
    integer :: buff_size

    write(this%name,'(a)') trim(filename)

    this%lock_time_dim_for_locations = lock_time_dim_for_locations

    if(present(comm))then
      istat = nf90_create(path=trim(filename),&
                          cmode=NF90_NETCDF4,&
                          ncid=this%ncid,&
                          comm=comm%MPI_VAL,&
                          info=MPI_INFO_NULL%MPI_VAL)
      if (istat /= NF90_NOERR) call handle_ncerr(istat, 'create_result_file: Error creating '//filename)
      this%parallel_access = .true.
    else
        istat = nf90_create(path=trim(filename),&
                          cmode=NF90_NETCDF4,&
                          ncid=this%ncid)
      if (istat /= NF90_NOERR) call handle_ncerr(istat, 'create_result_file: Error creating '//filename)
      this%parallel_access = .false.
    end if



    call add_global_meta_data(this%ncid)

    if (this%lock_time_dim_for_locations) then
      if(dim_t>0)then
        call init_temporal_dims(this%ncid, dim_t)
      else
        call init_temporal_dims(this%ncid)
      end if


      istat = nf90_inq_dimid(this%ncid, "n", dim_id_ulim)

      istat = nf90_def_var(this%ncid,name="step",xtype=NF90_SHORT,dimids=dim_id_ulim, varid=var_id_info)
      istat = nf90_put_att(this%ncid, var_id_info, 'description',  'initial=-1, advance=0, '&
                          //'forecast=1, analysis=2, unconstrained_analysis=3')
    end if

    istat = NF90_DEF_DIM(this%ncid, 'ensemble', n_modeltasks, dim_id_ens)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error defining dim ensemble')

    ! configuration file
    call GET_COMMAND_ARGUMENT(1,arg)
    istat = nf90_put_att(this%ncid, NF90_GLOBAL, 'TIE-GCM nml file path', trim(arg))
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error writing TIE-GCM nml file path')

    call read_text_file_into_buffer(arg,file_content_buff,buff_size)
    istat = nf90_put_att(this%ncid, NF90_GLOBAL, 'TIE-GCM nml file content', file_content_buff(1:buff_size))

    call GET_COMMAND_ARGUMENT(2,arg)
    istat = nf90_put_att(this%ncid, NF90_GLOBAL, 'PDAF nml file path', trim(arg))
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error writing PDAF nml file path')

    call read_text_file_into_buffer(arg,file_content_buff,buff_size)
    istat = nf90_put_att(this%ncid, NF90_GLOBAL, 'PDAF nml file content', file_content_buff(1:buff_size))

    ! success flag
    istat = nf90_def_var(this%ncid,name="exit_status",xtype=NF90_SHORT,varid=var_id_status)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error defining exit_status variable')
    istat = nf90_put_var(this%ncid, varid=var_id_status, values=1)
    istat = nf90_put_att(this%ncid, var_id_status, name='description', &
      values='zero if the program was terminated normally '&
      '(reached finalize routine). Else the program terminated '&
      'before the finalize routine due to some error or manual abort')

    ! run time
    istat = nf90_def_var(this%ncid,name="run_time",xtype=NF90_FLOAT,varid=var_id_runtime)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error defining run_time variable')
    istat = nf90_put_var(this%ncid, varid=var_id_runtime, values=-1)
    istat = nf90_put_att(this%ncid, var_id_runtime, name='long_name', &
      values='run time of this run')
    istat = nf90_put_att(this%ncid, var_id_runtime, name='units', &
      values='seconds')


  end subroutine

  !> Closes the file's NetCDF handle.
  subroutine result_file_close(this)

    implicit none

    ! arguments
    class(result_file) :: this

    ! local
    integer :: istat

    if(this%ncid>=0) then
      istat = nf90_close(this%ncid)
    end if

  end subroutine

  !> Writes the exit-status flag to this file and syncs it.
  subroutine result_file_set_exit_status(this,status)

    implicit none
    ! arguments
    class(result_file) :: this
    integer, intent(in) :: status

    ! local
    integer :: istat
    integer :: var_id_status

    if (this%ncid >= 0) then
      istat = nf90_inq_varid(this%ncid, 'exit_status', var_id_status)
      istat = nf90_put_var(this%ncid, varid=var_id_status, values=status)
    end if

    call this%sync()

  end subroutine

  !> Writes the accumulated run time to this file and syncs it.
  subroutine result_file_write_run_time(this)

    use mpi_module,only: time_totalrun

    implicit none
    ! arguments
    class(result_file) :: this

    ! local
    integer :: istat
    integer :: var_id_runtime

    if (this%ncid >= 0) then
      istat = nf90_inq_varid(this%ncid, 'run_time', var_id_runtime)
      istat = nf90_put_var(this%ncid, varid=var_id_runtime, values=time_totalrun)
    end if

    call this%sync()

  end subroutine

  !> Delegates to next_step for this file's root group, if open.
  subroutine result_file_next_step(this, algorithm_step, counter)

    implicit none

    ! arguments
    class(result_file) :: this
    character(len=*), intent(in) :: algorithm_step
    integer, intent(in) :: counter

    if(this%ncid>=0) then
      call next_step(this%ncid, algorithm_step, counter)
    end if

  end subroutine

  !> Flushes this file to disk.
  subroutine result_file_sync(this)

    implicit none

    ! arguments
    class(result_file) :: this

    ! local
    integer :: istat

    if(this%ncid>=0) then
      istat=nf90_sync(this%ncid)
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'result_file_writer_sync')
    end if

  end subroutine


  !> Initializes a spatial domain's name, member/moment flags, and write cadence.
  subroutine spatial_domain_init(this,domain_name,save_members,kmax,write_every_sec,&
                                 force_write_on_update,save_n_steps_after_update)

    ! tiegcm
    use cons_module,only: dt

    implicit none

    ! arguments
    class(spatial_domain), intent(inout) :: this
    character(len=*), intent(in) :: domain_name
    logical, intent(in) :: save_members
    integer, intent(in) :: kmax
    integer, intent(in) :: write_every_sec
    logical, intent(in) :: force_write_on_update
    integer, intent(in) :: save_n_steps_after_update

    this%domain_name = domain_name
    this%save_members = save_members
    this%save_moments = kmax>0
    this%kmax = kmax
    if(write_every_sec>0) then
      this%write_every = ceiling(real(write_every_sec) / dt)
    end if
    this%force_write_on_update = force_write_on_update
    this%save_n_steps_after_update = save_n_steps_after_update
  end subroutine

  !> Creates (if needed) the domain's NetCDF group in a file and assigns the requested mean/std/skew/kurt/members subgroup.
  subroutine spatial_domain_assign(this,file,group)
    ! assigns the given group (mean,std,members) at this spatial domain to a file

    ! interal
     use netcdf_functionality, only: init_temporal_dims

    implicit none

    ! arguments
    class(spatial_domain) :: this
    type(result_file), intent(in) :: file
    character(len=*), intent(in) :: group

    ! local
    integer :: istat
    integer :: ncid

    integer :: var_id_info, dim_id_ulim

    if (file%ncid >= 0) then

      istat = nf90_inq_ncid(file%ncid,this%domain_name,ncid)
      if (istat /= NF90_NOERR) then

        write(*,*) 'defining ', trim(this%domain_name), ' in ', file%name

        istat = nf90_def_grp(file%ncid,trim(this%domain_name),ncid)
        if (istat /= NF90_NOERR) call handle_ncerr(istat, 'spatial_domain_assign: Error creating '//group)

        if (file%lock_time_dim_for_locations.eqv..false.) then
          call init_temporal_dims(ncid)

          istat = nf90_inq_dimid(ncid, "n", dim_id_ulim)

          istat = nf90_def_var(ncid,name="step",xtype=NF90_SHORT,dimids=dim_id_ulim, varid=var_id_info)
          istat = nf90_put_att(ncid, var_id_info, 'description',  'initial=-1, advance=0, '&
                              //'forecast=1, analysis=2, unconstrained_analysis=3')
        end if

        call this%add_meta_data(ncid)
        call this%add_dimensions(ncid)

      end if
      select case(trim(group))
        case('members')
          write(*,*) 'defining ', trim(this%domain_name),'/members ', 'in ', file%name
          this%root_id_members = file%ncid
          this%spatial_id_members = ncid
          istat = nf90_def_grp(ncid, trim(group), this%grp_id_members)
          if (istat /= NF90_NOERR) call handle_ncerr(istat, 'spatial_domain_assign: Error creating members')
        case('mean')
          call spatial_domain_assign_moment(this,file,ncid,1,"ensemble mean")
        case('std')
          call spatial_domain_assign_moment(this,file,ncid,2,"unbiased standard deviation of ensemble")
        case('skew')
          call spatial_domain_assign_moment(this,file,ncid,3,"unbiased skewness of ensemble")
        case('kurt')
          call spatial_domain_assign_moment(this,file,ncid,4,"unbiased excess kurtosis of ensemble")
        case default
          ncid = -1 ! assign invalid id
      end select

    end if

  end subroutine

  !> Creates the NetCDF subgroup for one ensemble moment (mean/std/skew/kurt) and records its ids.
  subroutine spatial_domain_assign_moment(this, file, spatial_ncid, k, description)
    implicit none

    ! arguments
    class(spatial_domain) :: this
    type(result_file), intent(in) :: file
    integer, intent(in) :: spatial_ncid
    integer, intent(in) :: k
    character(len=*), intent(in) :: description

    ! local
    integer :: istat

    if(spatial_ncid>=0) then
      if(this%kmax>k-1) then
        write(*,*) 'defining ', trim(this%domain_name),' ',trim(moment_nc_name(k)), ' in ', file%name
        this%root_id_moments(k) = file%ncid
        this%spatial_id_moments(k) = spatial_ncid
        istat = nf90_def_grp(spatial_ncid, trim(moment_nc_name(k)), this%grp_id_moments(k))
        if (istat /= NF90_NOERR) call handle_ncerr(istat, 'spatial_domain_assign: Error creating '//trim(moment_nc_name(k)))
        istat = nf90_put_att(this%grp_id_moments(k), NF90_GLOBAL,"description",description)
      end if
    end if
  end subroutine

  !> Returns the NetCDF group id for a given data group name (members/mean/std/skew/kurt).
  function spatial_domain_get_grp_id(this,data_group) result(ncid)
    implicit none

    ! arguments
    class(spatial_domain) :: this
    character(len=*), intent(in) :: data_group
    ! returns
    integer :: ncid

    select case(data_group)
      case('members')
        ncid = this%grp_id_members
      case('mean')
        ncid = this%grp_id_moments(1)
      case('std')
        ncid = this%grp_id_moments(2)
      case('skew')
        ncid = this%grp_id_moments(3)
      case('kurt')
        ncid = this%grp_id_moments(4)
      case default
        ncid = -1 ! assign invalid id
    end select

  end function

  !> Returns the distinct root NetCDF ids among this domain's members/moment groups.
  subroutine spatial_domain_get_unique_spatial_ids (this, unique_roots)

    use m_unirnk, only: unirnk

    implicit none

    ! arguments
    class(spatial_domain) :: this
    integer, dimension(:), allocatable, intent(out) :: unique_roots

    ! local
    integer, dimension(5) :: ids
    integer, dimension(5) :: ordered
    integer :: n

    ids(1) = this%spatial_id_members
    ids(2:5) = this%spatial_id_moments

    call unirnk(ids,ordered,n)

    allocate(unique_roots(n))
    unique_roots = ids(ordered(1:n))
  end subroutine

  !> Advances the domain-local temporal counter/step flag when files don't share a global time dimension.
  subroutine spatial_domain_next_step(this, algorithm_step)

    ! tie-gcm
    use init_module, only: istep

    implicit none

    ! arguments
    class(spatial_domain) :: this
    character(len=*), intent(in) :: algorithm_step

    ! local
    integer, dimension(:), allocatable :: roots
    integer :: i

    if(this%result_writer%lock_time_dim_for_locations.eqv..false.) then
      call this%get_unique_spatial_ids(roots)

      this%counter = this%counter+1

      do i =1, size(roots)
        if(roots(i)>=0) then
          call next_step(roots(i), algorithm_step, this%counter)
        end if
      end do

      if (allocated(roots)) deallocate(roots)
    end if

    this%istep_last_write=istep

  end subroutine

  !> Returns the current time-step counter, from the writer or the domain itself depending on the time-dimension mode.
  function spatial_domain_get_counter(this) result(counter)

    implicit none

    ! arguments
    class(spatial_domain) :: this

    ! result
    integer :: counter

    if(this%result_writer%lock_time_dim_for_locations) then
      counter = this%result_writer%counter
    else
      counter = this%counter
    end if

  end function

  !> Decides whether the current step should be written, based on write cadence and algorithm step.
  function spatial_domain_time_to_write(this,algorithm_step) result(time_to_write)

    ! tie-gcm
    use init_module, only: istep

    implicit none

    ! arguments
    class(spatial_domain) :: this
    character(len=*), intent(in) :: algorithm_step

    ! result
    logical :: time_to_write

    ! local
    logical :: is_update_step
    logical :: is_within_range_after_update_step

    is_update_step =  (algorithm_step/="advance") .and. (algorithm_step/="initial")

    if(is_update_step) this%istep_last_update = istep

    is_within_range_after_update_step = .false.
    if(this%save_n_steps_after_update > 0) then
      if (this%istep_last_update >= 0) then
        if ((istep-this%istep_last_update >= 1) .and. &
           (istep-this%istep_last_update <= this%save_n_steps_after_update)) then
!             write(*,*) "Save step", istep, " since it is within the first ", this%save_n_steps_after_update, &
!                        " steps after the previous update at step", this%istep_last_update, &
!                        " You can control this via the output%save_n_steps_after_update option"
            is_within_range_after_update_step = .true.
        end if
      end if
    end if


    time_to_write = (modulo(istep,this%write_every)==0) .or. &                                       ! write every i-th step
                    (algorithm_step == "initial") .or. &                                             ! always write initial step
                    ((this%force_write_on_update.eqv..true.).and.(is_update_step.eqv..true.)) .or. & ! write all update steps (optional)
                    is_within_range_after_update_step                                                ! write n steps after update (optional)

    ! analysis step coincidents with this advance step. No need to write dublicate
    if((algorithm_step=="advance") .and. (istep==this%istep_last_write))then
      time_to_write = .false.
    end if

  end function

  !> Adds TIE-GCM spatial dimensions to a NetCDF group.
  subroutine nc_tiegcm_add_dimensions(this, ncid)

    ! intern
    use netcdf_functionality, only:  init_spacial_dims

    implicit none

    ! arguments
    class(nc_tiegcm), intent(inout) :: this
    integer, intent(in) :: ncid

    if ( ncid >= 0) then
      call init_spacial_dims(ncid)
    end if
  end subroutine

  !> Adds TIE-GCM-grid-specific metadata attributes (currently a no-op placeholder).
  subroutine nc_tiegcm_add_meta_data(this, ncid)

    implicit none

    ! arguments
    class(nc_tiegcm), intent(inout) :: this
    integer, intent(in) :: ncid

    if ( ncid >= 0) then

    end if
  end subroutine

  !> Defines a TIE-GCM-grid variable in the members and/or moment NetCDF groups.
  subroutine nc_tiegcm_add_variables(this, varname)

    ! intern
    use quantity_info_module, only: quantity_info, get_info

    implicit none

    ! arguments
    class(nc_tiegcm) :: this
    character(len=*), intent(in) :: varname

    ! local
    integer, dimension(5) :: var_dimensions
    integer :: var_id
    integer :: istat

    type(quantity_info), pointer :: info

    integer :: ncid

    integer :: i

    ! black list
    if(varname == 'P') then
      return
    end if

    info=>get_info(varname)

    if(this%save_members) then

      istat = nf90_inq_ncid(this%root_id_members,this%domain_name, ncid)
      if(istat==NF90_NOERR)then

        call nc_tiegcm_get_dimension_ids(ncid, info, var_dimensions, this%result_writer%lock_time_dim_for_locations)

        call rfw_def_var(ncid=this%grp_id_members,&
                        varname=varname,&
                        dimids=var_dimensions,&
                        var_id=var_id)

        call add_metadata_to_var(this%grp_id_members,var_id,varname)
      end if
    end if
    if(this%save_moments) then

      do i=1,this%kmax
        istat = nf90_inq_ncid(this%root_id_moments(i), this%domain_name, ncid)
        if(istat==NF90_NOERR)then
          call nc_tiegcm_get_dimension_ids( ncid, info, var_dimensions, this%result_writer%lock_time_dim_for_locations)

          call rfw_def_var(ncid=this%grp_id_moments(i),&
                          varname=varname,&
                          dimids=var_dimensions((/1,2,3,5/)),&
                          var_id=var_id)
          call add_metadata_to_var(this%grp_id_moments(i),var_id,varname)
        end if
      end do

    end if

  end subroutine

  !> Looks up the lon/lat/lev/ensemble/time dimension ids for a TIE-GCM NetCDF group.
  subroutine nc_tiegcm_get_dimension_ids(ncid, info, dimensions, lock_time_dim_for_locations)

    ! intern
    use quantity_info_module, only: quantity_info, LEVEL_INT,LEVEL_MID,LEVEL_NONE

    implicit none

    ! arguments
    integer, intent(in) :: ncid
    type(quantity_info), pointer, intent(in) :: info
    integer, dimension(5), intent(out) :: dimensions
    logical, intent(in) :: lock_time_dim_for_locations

    ! local
    integer :: parent_ncid
    integer :: istat
    integer :: ncid_n

    ! same order as in TIE-GCM history files
    istat = nf90_inq_dimid(ncid, "lon", dimensions(1))
    istat = nf90_inq_dimid(ncid, "lat", dimensions(2))
    select case(info%level)
        case(LEVEL_MID)
          istat = nf90_inq_dimid(ncid, "lev", dimensions(3))
        case(LEVEL_INT)
          istat = nf90_inq_dimid(ncid, "ilev", dimensions(3))
        case(LEVEL_NONE)
          ! degenerate vertical dimension of size one, so that quantities
          ! without vertical extent have the same rank as all others
          istat = nf90_inq_dimid(ncid, "lev1", dimensions(3))
        case default
          call shutdown('nc_tiegcm_get_dimension_ids: unhandled level of quantity '//trim(info%name))
    end select

    istat = nf90_inq_grp_parent(ncid,parent_ncid)

    istat = nf90_inq_dimid(parent_ncid, "ensemble", dimensions(4))

    if(lock_time_dim_for_locations) then
      ncid_n = parent_ncid
    else
      ncid_n = ncid
    end if

    istat = nf90_inq_dimid(ncid_n, "n", dimensions(5))

  end subroutine

  !> Reshapes and writes one TIE-GCM-grid field (member or moment) to NetCDF for the current step.
  subroutine nc_tiegcm_write(this,data_group,varname,data_tgcm_intern,data_state)

    ! tie-gcm
    use fields_module,only: levd0,lond0,latd0
    use mpi_module, only: mytid
    use params_module, only: nlon, nlat

    ! intern
    use mod_parallel_pdaf, only: task_id, filterpe
    use state_module,only: nlonX, nlatX, nlevX, idx_intern, idx_nc

    implicit none

    ! arguments
    class(nc_tiegcm) :: this
    character(len=*), intent(in) :: data_group
    ! assumed shape in the vertical, since quantities without vertical extent
    ! (LEVEL_NONE) are stored with a single level only
    real, dimension(levd0:,lond0:,latd0:), intent(in), optional :: data_tgcm_intern
    real, dimension(nlevX,nlonX,nlatX), intent(in), optional  :: data_state
    character(len=*), intent(in) :: varname

    ! local
    integer :: ncid
    integer :: var_id, istat
    integer :: count_p(5), start_p(5)
    real, dimension(:,:,:), allocatable :: field_reshaped

    if( (data_group/='members') .and. (filterpe.eqv..false.) ) then
      return
    end if

    ! black list
    if(varname == 'P') then
      return
    end if

    ! Order in netcdf is different then interal storage
    ! tiegcm: lev lon lat
    ! netcdf: lon lat lev

    count_p(1) = idx_nc(mytid)%nlons
    count_p(2) = idx_nc(mytid)%nlats


    if(present(data_tgcm_intern)) then
      ! one for quantities without vertical extent, nlevp1 otherwise
      count_p(3) = size(data_tgcm_intern,dim=1)
      allocate(field_reshaped(count_p(1),count_p(2),count_p(3)))
      field_reshaped(:,:,:) = reshape( data_tgcm_intern(:,idx_intern(mytid)%lon0:idx_intern(mytid)%lon1,&
                                                  idx_intern(mytid)%lat0:idx_intern(mytid)%lat1),&
                                count_p(1:3), order=(/3,1,2/))
    else if(present(data_state)) then
      count_p(3) = size(data_state,dim=1)
      allocate(field_reshaped(count_p(1),count_p(2),count_p(3)))
      field_reshaped(:,:,:) = reshape( data_state, &
                                count_p(1:3), order=(/3,1,2/))
    else
      write(*,*) "ERROR invalid call to nc_tiegcm_write"
    end if

    ncid = this%get_grp_id(data_group)

    select case(this%result_writer%strategy)
      case(SINGLE_FILE_SINGLE_WRITER)

        select case(data_group)
          case('members')
            call write_members_3d_lon_lat_distributed(&
              this,&
              ncid,&
              varname,&
              field_reshaped,&
              total_size=(/nlon,nlat,size(field_reshaped,dim=3)/),&
              offsets=(/idx_nc(mytid)%lon0,idx_nc(mytid)%lat0,levd0/))
          case('mean','std','skew','kurt')
            call write_moments_3d_lon_lat_distributed(&
              this,&
              ncid,&
              varname,&
              field_reshaped,&
              total_size=(/nlon,nlat,size(field_reshaped,dim=3)/),&
              offsets=(/idx_nc(mytid)%lon0,idx_nc(mytid)%lat0,levd0/))
        end select

      case(MULTIPLE_FILE_MULTIPLE_WRITER)

        istat = nf90_inq_varid(ncid, varname, var_id)
        if (istat /= NF90_NOERR) call handle_ncerr(istat,'nc_tiegcm_write: inq '//varname)

        start_p(1) = idx_nc(mytid)%lon0
        start_p(2) = idx_nc(mytid)%lat0
        start_p(3) = levd0

        select case(data_group)
          case('members')
            count_p(4) = 1
            count_p(5) = 1

            start_p(4) = task_id

            start_p(5) = this%get_counter()


            istat = nf90_put_var(ncid, var_id, start=start_p, count=count_p, values=field_reshaped)
            if (istat /= NF90_NOERR) call handle_ncerr(istat,'nc_tiegcm_write: put members of '//varname)

          case('mean','std','skew','kurt')
            if(filterpe) then
              count_p(4) = 1
              start_p(4) = this%get_counter()

              istat = nf90_put_var(ncid, var_id, start=start_p(1:4), count=count_p(1:4), values=field_reshaped)
              if (istat /= NF90_NOERR) call handle_ncerr(istat,'nc_tiegcm_write: put moments of '//varname)

            end if
        end select

    end select

    if(allocated(field_reshaped)) deallocate(field_reshaped)
  end subroutine

  !> Adds the regular-grid dataset's spatial dimensions to a NetCDF group.
  subroutine nc_reg_grid_add_dimensions(this, ncid)

    implicit none
    ! arguments
    class(nc_reg_grid), intent(inout) :: this
    integer, intent(in) :: ncid

    if ( ncid >= 0) then
       call this%dataset%add_spatial_dim(ncid)
    end if

  end subroutine

  !> Copies grid-file provenance metadata (creation date, source path) into the output group.
  subroutine nc_reg_grid_add_meta_data(this, ncid)

    implicit none
    ! arguments
    class(nc_reg_grid), intent(inout) :: this
    integer, intent(in) :: ncid

    ! local
    integer :: istat
    character(len=32) :: date_created

    if ( ncid >= 0) then

      istat = nf90_get_att(this%dataset%ncid, NF90_GLOBAL, 'date_created', date_created)
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error getting date_created attribute')

      istat = nf90_put_att(ncid, &
        NF90_GLOBAL,"grid_file_date_created", &
        trim(date_created))
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error writing attribute grid_file_date_created')

      istat = nf90_put_att(ncid, &
        NF90_GLOBAL,"grid_file_path", &
        trim( this%dataset%nc_file))
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error writing attribute grid_file_path')
    end if
  end subroutine

  !> Defines a regular-grid variable in the members and/or moment NetCDF groups.
  subroutine nc_reg_grid_add_variables(this, varname)

    implicit none
    ! arguments
    class(nc_reg_grid) :: this
    character(len=*), intent(in) :: varname

    ! local
    integer, dimension(5) :: var_dimensions
    integer :: var_id
    integer :: istat

    integer :: ncid

    integer :: i

    if(this%save_members) then

      istat = nf90_inq_ncid(this%root_id_members, this%domain_name, ncid)
      if(istat==NF90_NOERR)then
        call nc_reg_grid_get_dimension_ids(ncid, var_dimensions, this%result_writer%lock_time_dim_for_locations)

        call rfw_def_var(ncid=this%grp_id_members,&
                          varname=varname,&
                          dimids=var_dimensions,&
                          var_id=var_id)
        call add_metadata_to_var(this%grp_id_members,var_id,varname)
      end if
    end if

    if(this%save_moments) then

      do i=1,this%kmax
        istat = nf90_inq_ncid(this%root_id_moments(i), this%domain_name, ncid)
        if(istat==NF90_NOERR)then
          call nc_reg_grid_get_dimension_ids( ncid, var_dimensions, this%result_writer%lock_time_dim_for_locations)

          call rfw_def_var(ncid=this%grp_id_moments(i),&
                            varname=varname,&
                            dimids=var_dimensions((/1,2,3,5/)),&
                            var_id=var_id)
          call add_metadata_to_var(this%grp_id_moments(i),var_id,varname)
        end if
      end do

    end if

  end subroutine

  !> Looks up the lon/lat/alt/ensemble/time dimension ids for a regular-grid NetCDF group.
  subroutine nc_reg_grid_get_dimension_ids(ncid, dimensions,lock_time_dim_for_locations)

  implicit none

    ! arguments
    integer, intent(in) :: ncid
    integer, dimension(5), intent(out) :: dimensions
    logical, intent(in) :: lock_time_dim_for_locations

    ! local
    integer :: parent_ncid
    integer :: istat
    integer :: ncid_n

    istat = nf90_inq_dimid(ncid, "lon", dimensions(1))
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'nc_reg_grid_add_variables inquire dimension lon')
    istat = nf90_inq_dimid(ncid, "lat", dimensions(2))
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'nc_reg_grid_add_variables inquire dimension lat')
    istat = nf90_inq_dimid(ncid, "alt", dimensions(3))
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'nc_reg_grid_add_variables inquire dimension alt')

    istat = nf90_inq_grp_parent(ncid,parent_ncid)

    istat = nf90_inq_dimid(parent_ncid, "ensemble", dimensions(4))
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'nc_reg_grid_add_variables inquire dimension ensemble')

    if(lock_time_dim_for_locations) then
      ncid_n = parent_ncid
    else
      ncid_n = ncid
    end if

    istat = nf90_inq_dimid(ncid_n, "n", dimensions(5))
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'nc_reg_grid_add_variables inquire dimension time')

  end subroutine

  !> Writes one regular-grid field (member or moment) to NetCDF for the current step.
  subroutine nc_reg_grid_write(this,data_group,varname,data)

    ! intern
    use mod_parallel_pdaf, only: task_id, filterpe

    implicit none

    ! arguments
    class(nc_reg_grid) :: this
    character(len=*), intent(in) :: data_group
    real, dimension(:,:,:), intent(in)  :: data
    character(len=*), intent(in) :: varname

    ! local
    integer :: ncid
    integer nlon, nlat, nalt

    integer, dimension(5) :: start_p, count_p
    integer :: istat
    integer :: var_id


    if( (data_group/='members') .and. (filterpe.eqv..false.) ) then
      return
    end if

    ncid = this%get_grp_id(data_group)

    nlon=size(this%dataset%lon)
    nlat=size(this%dataset%lat)
    nalt = (this%dataset%alt_last-this%dataset%alt_first)+1

    select case(this%result_writer%strategy)
      case(SINGLE_FILE_SINGLE_WRITER)

        select case(data_group)
          case('members')
            call write_members_3d_lon_lat_distributed(&
              this,&
              ncid,&
              varname,&
              data,&
              total_size=(/nlon, nlat, nalt/),&
              offsets=(/this%dataset%offset_lon,this%dataset%offset_lat,1/),&
              vert_start=this%dataset%alt_first)
          case('mean','std','skew','kurt')
            call write_moments_3d_lon_lat_distributed(&
              this,&
              ncid,&
              varname,&
              data,&
              total_size=(/nlon, nlat, nalt/),&
              offsets=(/this%dataset%offset_lon,this%dataset%offset_lat,1/),&
              vert_start=this%dataset%alt_first)
        end select
      case(MULTIPLE_FILE_MULTIPLE_WRITER)
        istat = nf90_inq_varid(ncid, varname, var_id)
        if (istat /= NF90_NOERR) call handle_ncerr(istat,'nc_reg_grid_write: inq '//varname)

        start_p(1) = this%dataset%offset_lon
        start_p(2) = this%dataset%offset_lat
        start_p(3) = this%dataset%alt_first

        count_p(1) = this%dataset%nlon
        count_p(2) = this%dataset%nlat
        count_p(3) = nalt

        select case(data_group)
          case('members')
            count_p(4) = 1
            count_p(5) = 1

            start_p(4) = task_id
            start_p(5) = this%get_counter()

            istat = nf90_put_var(ncid, var_id, start=start_p, count=count_p, values=data)
          case('mean','std','skew','kurt')

            count_p(4) = 1
            start_p(4) = this%get_counter()

            if(filterpe) then
              istat = nf90_put_var(ncid, var_id, start=start_p(1:4), count=count_p(1:4), values=data)
            end if
        end select
        if (istat /= NF90_NOERR) call handle_ncerr(istat,'nc_reg_grid_write: put '//varname)
    end select
  end subroutine

  !> Associates a regular-grid dataset with this spatial domain.
  subroutine nc_reg_grid_link_dataset(this, dataset)

    implicit none

    ! arguments
    class(nc_reg_grid) :: this
    type(reg_grid_dataset_root), target :: dataset

    this%dataset => dataset
  end subroutine

  !> Defines the lon/lat/alt position variables for a trajectory NetCDF group.
  subroutine nc_trajectory_add_dimensions(this, ncid)

    implicit none

    ! arguments
    class(nc_trajectory), intent(inout) :: this
    integer, intent(in) :: ncid

    ! local
    integer :: parent_ncid
    integer :: dim_id(2)
    integer :: istat
    integer :: var_id

    if ( ncid >= 0) then

      call nc_trajectory_get_dimension_ids(ncid, dim_id, this%result_writer%lock_time_dim_for_locations)

      call rfw_def_var(ncid,'lon',dim_id(2:2),var_id)
      istat = nf90_put_att(ncid,var_id,"long_name", "geographic longitude (-west, +east)")
      istat = nf90_put_att(ncid,var_id,"units",'degrees_east')

      call rfw_def_var(ncid,'lat',dim_id(2:2),var_id)
      istat = nf90_put_att(ncid,var_id,"long_name", "geocentric latitude (-south, +north)")
      istat = nf90_put_att(ncid,var_id,"units",'degrees_east')

      call rfw_def_var(ncid,'alt',dim_id(2:2),var_id)
      istat = nf90_put_att(ncid,var_id,"long_name", "altitude")
      istat = nf90_put_att(ncid,var_id,"units",'m')
    end if
  end subroutine

  !> Sets the CF featureType="trajectory" attribute on the group.
  subroutine nc_trajectory_add_meta_data(this, ncid)

    implicit none

    ! arguments
    class(nc_trajectory), intent(inout) :: this
    integer, intent(in) :: ncid

    ! local
    integer :: istat

    if ( ncid >= 0) then
       istat = nf90_put_att(ncid, NF90_GLOBAL,"featureType", "trajectory")
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error writing attribute featureType')
    end if
  end subroutine

  !> Defines a scalar trajectory variable in the members and/or moment NetCDF groups.
  subroutine nc_trajectory_add_variables(this, varname)

    ! intern
    use quantity_info_module, only: supports_point_interpolation

    implicit none

    ! arguments
    class(nc_trajectory) :: this
    character(len=*), intent(in) :: varname

    ! local
    integer, dimension(2) :: var_dimensions
    integer :: var_id
    integer :: istat

    integer :: ncid

    integer :: i

    ! quantities that cannot be interpolated to a point are not part of a
    ! trajectory file
    if(supports_point_interpolation(varname).eqv..false.) then
      write(*,*) 'INFO: ', trim(varname), ' cannot be interpolated to a point and is', &
                 ' therefore not written to trajectory ', trim(this%domain_name)
      return
    end if

    if(this%save_members) then

      istat = nf90_inq_ncid(this%root_id_members, this%domain_name, ncid)
      if(istat==NF90_NOERR)then
        call nc_trajectory_get_dimension_ids( ncid, var_dimensions, this%result_writer%lock_time_dim_for_locations)

        call rfw_def_var(ncid=this%grp_id_members,&
                          varname=varname,&
                          dimids=var_dimensions,&
                          var_id=var_id)
        call add_metadata_to_var(this%grp_id_members,var_id,varname)
      end if
    end if

    if(this%save_moments) then

      do i=1,this%kmax
        istat = nf90_inq_ncid(this%root_id_moments(i), this%domain_name, ncid)
        if(istat==NF90_NOERR)then
          call nc_trajectory_get_dimension_ids( ncid, var_dimensions, this%result_writer%lock_time_dim_for_locations)

          call rfw_def_var(ncid=this%grp_id_moments(i),&
                            varname=varname,&
                            dimids=var_dimensions(2:2),&
                            var_id=var_id)
          call add_metadata_to_var(this%grp_id_moments(i),var_id,varname)
        end if
      end do

    end if

  end subroutine

  !> Looks up the ensemble/time dimension ids for a trajectory NetCDF group.
  subroutine nc_trajectory_get_dimension_ids(ncid, dimensions,lock_time_dim_for_locations)

  implicit none

    ! arguments
    integer, intent(in) :: ncid
    integer, dimension(2), intent(out) :: dimensions
    logical, intent(in) :: lock_time_dim_for_locations

    ! local
    integer :: parent_ncid
    integer :: istat

    integer :: ncid_n

    istat = nf90_inq_grp_parent(ncid,parent_ncid)

    istat = nf90_inq_dimid(parent_ncid, "ensemble", dimensions(1))
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'nc_trajectory_add_variables inquire dimension ensemble')

    if(lock_time_dim_for_locations) then
      ncid_n = parent_ncid
    else
      ncid_n = ncid
    end if
    istat = nf90_inq_dimid(ncid_n, "n", dimensions(2))
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'nc_trajectory_add_variables inquire dimension time')

  end subroutine

  !> Gathers and writes one scalar trajectory value (member or moment) to NetCDF for the current step.
  subroutine nc_trajectory_write(this,data_group,varname,val)

    ! ATTENTION assume val is located on root
    ! nc_trajectory is not using parallel io

    ! tie-gcm
    use mpi_module, only: mytid

    ! intern
    use mod_parallel_pdaf, only: filterpe, n_modeltasks, rank_world, task_id

    implicit none

    ! arguments
    class(nc_trajectory) :: this
    character(len=*), intent(in) :: data_group
    real, intent(in) :: val
    character(len=*), intent(in) :: varname

    ! local
    integer :: ncid
    integer :: var_id
    integer :: istat

    real, dimension(1) :: val_copy
    real, dimension(:), allocatable :: members

    if( (data_group/='members') .and. (filterpe.eqv..false.) ) then
      return
    end if

    ncid = this%get_grp_id(data_group)
    istat = nf90_inq_varid(ncid, varname, var_id)

    val_copy = val

        select case(data_group)
          case('members')
!             call collect_point_on_root(on_this_rank,val_copy(1))

            if ( mytid == 0) then
              allocate(members(n_modeltasks))
            else
              allocate(members(0))
            end if

            if ( mytid == 0) then
              call gather_ensemble_1d(val_copy(1), members)
            end if
            if ( rank_world == 0) then
              ! ATTENTION unbuffered writing is slow
              istat = nf90_put_var(ncid,&
                                   var_id,&
                                   start=(/1,this%get_counter()/),&
                                   count=(/n_modeltasks,1/),&
                                   values=members)
            end if

            deallocate(members)
          case('mean','std','skew','kurt')
            if(filterpe) then
!                call collect_point_on_root(on_this_rank,val_copy(1))
               if ( mytid == 0) then

                istat = nf90_put_var(ncid=ncid,&
                          varid=var_id,&
                          values=val_copy,&
                          count=(/1/),&
                          start=(/this%get_counter()/))
                end if
              end if
        end select
  end subroutine

  !> Writes the current lon/lat/alt satellite position to all root NetCDF groups of this trajectory.
  subroutine nc_trajectory_write_position(this, position)

    ! intern
    use mod_parallel_pdaf, only: rank_world

    implicit none

    ! arguments
    class(nc_trajectory) :: this
    real, dimension(3), intent(in) :: position

    ! local
    integer :: istat
    integer, dimension(:), allocatable :: roots
    integer :: i
    integer :: var_id

    call this%get_unique_spatial_ids(roots)

    if ( rank_world == 0) then
      do i=1,size(roots)
        if(roots(i)>=0) then
          istat = nf90_inq_varid(roots(i), 'lon', var_id)
          if (istat /= NF90_NOERR) call handle_ncerr(istat,'nc_trajectory_write_position inquire lon')
          istat = nf90_put_var(ncid=roots(i),&
                               varid=var_id,&
                               values=position(1:1),&
                               count=(/1/),&
                               start=(/this%get_counter()/))
          if (istat /= NF90_NOERR) call handle_ncerr(istat,'nc_trajectory_write_position put lon')

          istat = nf90_inq_varid(roots(i), 'lat', var_id)
          if (istat /= NF90_NOERR) call handle_ncerr(istat,'nc_trajectory_write_position inquire lat')
          istat = nf90_put_var(ncid=roots(i),&
                               varid=var_id,&
                               values=position(2:2),&
                               count=(/1/),&
                               start=(/this%get_counter()/))
          if (istat /= NF90_NOERR) call handle_ncerr(istat,'nc_trajectory_write_position put lat')

          istat = nf90_inq_varid(roots(i), 'alt', var_id)
          if (istat /= NF90_NOERR) call handle_ncerr(istat,'nc_trajectory_write_position inquire alt')
          istat = nf90_put_var(ncid=roots(i),&
                               varid=var_id,&
                               values=position(3:3),&
                               count=(/1/),&
                               start=(/this%get_counter()/))
          if (istat /= NF90_NOERR) call handle_ncerr(istat,'nc_trajectory_write_position put alt')
        end if
      end do
    end if

    deallocate(roots)

  end subroutine

  !> Associates a trajectory dataset with this spatial domain.
  subroutine nc_trajectory_link_dataset(this, dataset)

    implicit none

    ! arguments
    class(nc_trajectory) :: this
    type(trajectory_data), target, intent(in) :: dataset

    this%dataset => dataset
  end subroutine

  !> Gathers a distributed 3D lon/lat field onto the root and writes it as one ensemble-moment slice.
  subroutine write_moments_3d_lon_lat_distributed(this,ncid,varname,data,total_size,offsets,vert_start)

    ! tie-gcm
    use mpi_module, only: mytid

    ! intern
    use mod_parallel_pdaf, only: filterpe

    implicit none

    ! arguments
    class(spatial_domain) :: this
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: varname
    real, dimension(:,:,:), intent(in) :: data
    integer, dimension(3), intent(in) :: total_size ! size of gatherd array only relevant at root
    integer, dimension(3), intent(in) :: offsets ! offsets of subdomains
    integer, optional :: vert_start ! start index of vertical dimension in netcdf file

    ! local
    real, dimension(:,:,:), allocatable :: member
    integer :: istat
    integer :: var_id

    integer :: start_p(4)
    integer :: count_p(4)

    ! filterpe == .true. and mytid == 0 ---> rank_world == 0
    if(filterpe) then
      if ( mytid == 0) then
        allocate(member(total_size(1), total_size(2), total_size(3)))
      else
        allocate(member(0,0,0))
      end if

      call gather_model_3d(data,&
                            member, &
                            (/offsets(1),offsets(2),offsets(3)/))
      if ( mytid == 0) then

        istat = nf90_inq_varid(ncid, varname, var_id)
        if (istat /= NF90_NOERR) call handle_ncerr(istat,'write_moments_3d_lon_lat_distributed: inq '//varname)

        start_p(1:3) = 1
        start_p(4) = this%get_counter()
        if( present(vert_start)) start_p(3)=vert_start

        count_p(1:3)=shape(member)
        count_p(4)=1

        istat = nf90_put_var(ncid, var_id, start=start_p, count=count_p, values=member)
        if (istat /= NF90_NOERR) call handle_ncerr(istat,'write_moments_3d_lon_lat_distributed: put '//varname)
      end if

      if (allocated(member)) deallocate(member)

    end if

  end subroutine

  !> Gathers a distributed 3D field across subdomains and ensemble members, then writes it as the members slice.
  subroutine write_members_3d_lon_lat_distributed(this,ncid,varname,data,total_size,offsets,vert_start)

    ! tie-gcm
    use mpi_module, only: mytid

    ! intern
    use mod_parallel_pdaf, only: rank_world, n_modeltasks

    implicit none

    ! arguments
    class(spatial_domain) :: this
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: varname
    real, dimension(:,:,:), intent(in) :: data

    integer, dimension(3), intent(in) :: total_size ! size of gatherd array only relevant at root
    integer, dimension(3), intent(in) :: offsets ! offsets of subdomains
    integer, optional :: vert_start ! start index of vertical dimension in netcdf file

    ! local
    real, dimension(:,:,:), allocatable :: member
    real, dimension(:,:,:,:), allocatable :: ensemble

    integer :: istat
    integer :: var_id

    integer :: start_p(5)
    integer :: count_p(5)


    if ( mytid == 0) then
      allocate(member(total_size(1), total_size(2), total_size(3)))
    else
      allocate(member(0,0,0))
    end if
    call gather_model_3d(data,&
                          member, &
                          (/offsets(1),offsets(2),offsets(3)/))

    if ( rank_world == 0) then
      allocate(ensemble(total_size(1), total_size(2), total_size(3), n_modeltasks))
    else
      allocate(ensemble(0,0,0,0))
    end if
    call gather_ensemble_3d(member, ensemble)

    if ( rank_world == 0) then

      istat = nf90_inq_varid(ncid, varname, var_id)
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'write_members_3d_lon_lat_distributed: inq '//varname)

      start_p(1:4) = 1
      start_p(5) = this%get_counter()
      if( present(vert_start)) start_p(3)=vert_start

      count_p(1:4)=shape(ensemble)
      count_p(5)=1

      istat = nf90_put_var(ncid, var_id, start=start_p, count=count_p, values=ensemble)
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'write_members_3d_lon_lat_distributed: put '//varname)
    end if

    if (allocated(member)) deallocate(member)
    if (allocated(ensemble)) deallocate(ensemble)

  end subroutine

  !> MPI-gathers each rank's 3D subdomain into a single global array on rank 0.
  subroutine gather_model_3d(send, recv, start_p)

    ! extern
    use mpi_f08

    ! tie-gcm
    use mpi_module, only: mytid, TIEGCM_WORLD, ntask, handle_mpi_err

    ! intern
    use array_print_module, only: printMat

    implicit none

    ! arguments
    real, dimension(:,:,:), intent(in) :: send ! subdomain of model
    real, dimension(:,:,:), intent(inout) :: recv

    integer, dimension(3), intent(in) :: start_p

    ! local
    integer :: ierr
    integer :: n

    integer :: starts(3,0:ntask-1)
    integer :: counts(3,0:ntask-1)

    integer :: upper_bound(3,0:ntask-1)
    integer :: tid

    character(len=128) :: err_buff

    ! collect start index and size of all ranks at root
    call MPI_Gather(sendbuf=start_p,&
                    sendcount=3,&
                    sendtype=MPI_Integer,&
                    recvbuf=starts,&
                    recvcount=3,&
                    recvtype=MPI_Integer,&
                    root=0,&
                    comm=TIEGCM_WORLD,&
                    ierror=ierr)
    call MPI_Gather(sendbuf=shape(send),&
                    sendcount=3,&
                    sendtype=MPI_Integer,&
                    recvbuf=counts,&
                    recvcount=3,&
                    recvtype=MPI_Integer,&
                    root=0,&
                    comm=TIEGCM_WORLD,&
                    ierror=ierr)


!     call printMat(send,"send",dim=1)
    if ( mytid == 0) then

      upper_bound = starts + counts - 1

      ! consitency check
      if(sum(product(counts,dim=1)) /= size(recv)) then
        write(*,*) "consitency check  failed"
        call printMat(starts,"lower bound")
        call printMat(upper_bound,"upper bound")
        call printMat(counts,"counts")
        write(err_buff,*) 'gather_model_3d: size of sends != size of receive', &
          sum(product(counts,dim=1)) ," /= ", size(recv)
        call shutdown(trim(err_buff))
      end if

      ! send data at root is written directly to recv
      recv(starts(1,mytid):upper_bound(1,mytid),&
           starts(2,mytid):upper_bound(2,mytid),&
           starts(3,mytid):upper_bound(3,mytid)) = send

      do tid=1,ntask-1
        n = product(counts(:,tid))

        call MPI_Recv(&
          buf=recv(starts(1,tid):upper_bound(1,tid),&
                   starts(2,tid):upper_bound(2,tid),&
                   starts(3,tid):upper_bound(3,tid)), &
          count=n, &
          datatype=MPI_REAL8, &
          source=tid,&
          tag=tid,&
          comm=TIEGCM_WORLD,&
          status=MPI_STATUS_IGNORE,&
          ierror=ierr)
        if (ierr /= 0) call handle_mpi_err(ierr,'MPI_Recv')
      end do

!       call printMat(recv,"gatherd",dim=1)
    else
      call MPI_Send(&
        buf=send, &
        count=size(send), &
        datatype=MPI_REAL8, &
        dest=0,&
        tag=mytid,&
        comm=TIEGCM_WORLD,&
        ierror=ierr)
      if (ierr /= 0) call handle_mpi_err(ierr,'MPI_Send')
    end if

  end subroutine

  !> MPI-gathers one 3D field from every ensemble task onto the coupling root.
  subroutine gather_ensemble_3d(send, recv)

    ! extern
    use mpi_f08

    ! tie-gcm
    use mpi_module, only: mytid

    ! intern
    use mod_parallel_pdaf, only: COMM_couple

    implicit none

    ! arguments
    real, dimension(:,:,:), intent(in) :: send ! model
    real, dimension(:,:,:,:), intent(inout) :: recv ! last dim is ensemble

    ! local
    integer :: ierr

    if(mytid==0) then
      call MPI_Gather(sendbuf=send,&
                  sendcount=size(send),&
                  sendtype=MPI_REAL8,&
                  recvbuf=recv,&
                  recvcount=size(send),&
                  recvtype=MPI_REAL8,&
                  root=0,&
                  comm=COMM_couple,&
                  ierror=ierr)
    end if

  end subroutine

  !> MPI-gathers one scalar value from every ensemble task onto the coupling root.
  subroutine gather_ensemble_1d(send, recv)

    ! extern
    use mpi_f08

    ! tie-gcm
    use mpi_module, only: mytid

    ! intern
    use mod_parallel_pdaf, only: COMM_couple

    implicit none

    ! arguments
    real, intent(in) :: send ! model
    real, dimension(:), intent(inout) :: recv ! last dim is ensemble

    ! local
    integer :: ierr

    if(mytid==0) then
      call MPI_Gather(sendbuf=send,&
                  sendcount=1,&
                  sendtype=MPI_REAL8,&
                  recvbuf=recv,&
                  recvcount=1,&
                  recvtype=MPI_REAL8,&
                  root=0,&
                  comm=COMM_couple,&
                  ierror=ierr)
    end if

  end subroutine

  !> Sends a scalar value from one designated rank to rank 0 via a point-to-point MPI message.
  subroutine collect_point_on_root(on_this_rank,val)

    ! extern
    use mpi_f08

    ! tie-gcm
    use mpi_module, only: mytid, TIEGCM_WORLD

    implicit none

    ! arguments
    logical, intent(in) :: on_this_rank
    real, intent(inout) :: val

    ! local
    integer :: ier
    TYPE(MPI_Status) :: status

    if(on_this_rank)then
     call mpi_send(val, 1, MPI_REAL8, &
                    0, 42, TIEGCM_WORLD, ier)
    end if

    if(mytid==0) then
      call mpi_recv(val, 1, MPI_REAL8, &
                     MPI_ANY_SOURCE, 42, TIEGCM_WORLD, status, ier)
    end if


  end subroutine

  !> Writes the units and long_name attributes for a NetCDF variable from the quantity registry.
  subroutine add_metadata_to_var(ncid,var_id,varname)

    ! intern
    use quantity_info_module, only: quantity_info, get_info

    implicit none
    ! arguments
    integer, intent(in) :: ncid
    integer, intent(in) :: var_id
    character(len=*), intent(in) :: varname

    ! local
    type(quantity_info), pointer :: info
    integer :: istat

    info=>get_info(varname)

    istat = nf90_put_att(ncid,var_id, "units", trim(adjustl(info%unit)))
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error putting units attribute')

    istat = nf90_put_att(ncid,var_id, "long_name", trim(adjustl(info%long_name)))
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error putting units attribute')

  end subroutine

  !> Creates the members/mean/std/skew/kurt NetCDF subgroups under a given group.
  subroutine add_member_and_moments_group(ncid,grp_id_members,grp_id_moments)

    implicit none

    ! arguments
    integer, intent(in) :: ncid
    integer, intent(out), optional :: grp_id_members
    integer, dimension(:), intent(out), optional :: grp_id_moments

    ! local
    integer :: istat
    integer :: kmax

    if(present(grp_id_members)) then
      istat = nf90_def_grp(ncid, 'members', grp_id_members)
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'add_member_and_moments_group members')
    end if
    if(present(grp_id_moments)) then
      kmax = size(grp_id_moments,dim=1)
      if(kmax>0)then
        istat = nf90_def_grp(ncid, 'mean', grp_id_moments(1))
        if (istat /= NF90_NOERR) call handle_ncerr(istat,'add_member_and_moments_group mean')
      end if
      if(kmax>1)then
        istat = nf90_def_grp(ncid, 'std', grp_id_moments(2))
        if (istat /= NF90_NOERR) call handle_ncerr(istat,'add_member_and_moments_group std')
      end if
      if(kmax>2)then
        istat = nf90_def_grp(ncid, 'skew', grp_id_moments(3))
        if (istat /= NF90_NOERR) call handle_ncerr(istat,'add_member_and_moments_group skew')
      end if
      if(kmax>3)then
        istat = nf90_def_grp(ncid, 'kurt', grp_id_moments(4))
        if (istat /= NF90_NOERR) call handle_ncerr(istat,'add_member_and_moments_group kurt')
      end if
    end if

  end subroutine

  !> Defines a compressed NetCDF variable with the module's standard type/chunking/deflate settings.
  subroutine rfw_def_var(ncid,varname,dimids,var_id)

    implicit none

    ! arguments
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: varname
    integer, dimension(:), intent(in) :: dimids
    integer, intent(out) :: var_id

    ! local
    integer :: istat
    character(len=64) :: err_buff

    istat = NF90_DEF_VAR(ncid=ncid,&
                             name=trim(varname),&
                             xtype=nc_var_xtype,&
                             dimids=dimids,&
                             varid=var_id,&
                             shuffle = nc_var_shuffle,&
                             fletcher32 = nc_var_fletcher32,&
                             deflate_level = nc_var_deflate_level )
    if (istat /= NF90_NOERR) then
      write(err_buff,'(a,a,a,i8)') 'rfw_def_var: ', varname, ' ncid:', ncid
      call handle_ncerr(istat, err_buff)
    end if

  end subroutine

  !> Reads an entire text file into a single character buffer (used to embed namelist files as NetCDF attributes).
  subroutine read_text_file_into_buffer(file_path,file_content_buffer, count)

    implicit none

    ! arguments
    character(len=*), intent(in) :: file_path
    character(len=*), intent(inout) :: file_content_buffer
    integer, intent(out) :: count

    ! local
    character(len=256) :: line_buff
    integer :: istat
    integer :: last, first
    integer, parameter :: fid = 53

    istat = 0
    first=1
    last = 0
    open(unit=fid,file=file_path,status='old')
    do while (istat == 0)
        read(fid,'(a)',iostat=istat) line_buff
        last = first + len_trim(line_buff) ! we do not need to substract 1 since linebreak char is added
        file_content_buffer(first:last) = trim(line_buff)//ACHAR(10)
        first=last+1
    end do
    close(fid)

    count = last

  end subroutine

  !> Immediately writes one TIE-GCM-grid field to an open NetCDF file/group, creating dims/vars on demand (for debugging).
  subroutine instant_tiegcm_write(ncid,name,level,t_idx,data_tgcm_intern,data_state)
    ! standalone subroutine to write data on tiegcm grid to nc file

    ! tie-gcm
    use fields_module,only: levd0,levd1,lond0,lond1,latd0,latd1
    use mpi_module, only: mytid
    use params_module, only: nlon, nlat, nlevp1

    ! intern
    use netcdf_functionality, only: init_temporal_dims, init_spacial_dims, add_model_time
    use quantity_info_module, only: LEVEL_INT,LEVEL_MID
    use state_module,only: nlonX, nlatX, nlevX, idx_intern, idx_nc

    implicit none

    ! arguments
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer, intent(in) :: level
    integer, intent(in) :: t_idx
    real, dimension(levd0:levd1,lond0:lond1,latd0:latd1), intent(in), optional :: data_tgcm_intern
    real, dimension(nlevX,nlonX,nlatX), intent(in), optional  :: data_state

    ! local
    integer :: istat
    integer :: var_dimensions(4)
    integer :: var_id
    integer :: len_ulim

    integer :: subdomain_size(3)
    real, dimension(:,:,:), allocatable :: field_reshaped
    real, dimension(:,:,:), allocatable :: member

    integer :: start_p(4)
    integer :: count_p(4)

    if ( mytid == 0) then

      ! create dimensions if missing
      istat = nf90_inq_dimid(ncid, "n", var_dimensions(4))
      if (istat /= NF90_NOERR) then
        call init_temporal_dims(ncid)
      end if

      istat = nf90_inq_dimid(ncid, "lon", var_dimensions(1))
      if (istat /= NF90_NOERR) then
        call init_spacial_dims(ncid)
      end if

      ! get dimension ids
      istat = nf90_inq_dimid(ncid, "lon", var_dimensions(1))
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'instant_tiegcm_write: inquire lon dimension')
      istat = nf90_inq_dimid(ncid, "lat", var_dimensions(2))
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'instant_tiegcm_write: inquire lat dimension')
      select case(level)
          case(LEVEL_MID)
            istat = nf90_inq_dimid(ncid, "lev", var_dimensions(3))
            if (istat /= NF90_NOERR) call handle_ncerr(istat,'instant_tiegcm_write: inquire lev dimension')
          case(LEVEL_INT)
            istat = nf90_inq_dimid(ncid, "ilev", var_dimensions(3))
            if (istat /= NF90_NOERR) call handle_ncerr(istat,'instant_tiegcm_write: inquire ilev dimension')
      end select

      istat = nf90_inq_dimid(ncid, "n", var_dimensions(4))
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'instant_tiegcm_write: inquire unlimited dimension')

      ! expand time dimension if necessary
      istat = nf90_inquire_dimension(ncid, var_dimensions(4), len=len_ulim)
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'instant_tiegcm_write: inquire length of dimension n')
      if (len_ulim < t_idx) then
        call add_model_time(ncid, t_idx, .false.)
      end if

      ! if variable is not existing create it
      istat = nf90_inq_varid(ncid, trim(name), var_id)
      if (istat /= NF90_NOERR) then
        call rfw_def_var(ncid=ncid,&
                         varname=name,&
                         dimids=var_dimensions,&
                         var_id=var_id)
      end if
    end if

    ! reshape the tiegcm to netcdf order
    if(present(data_tgcm_intern)) then
      subdomain_size = (/idx_nc(mytid)%nlons,idx_nc(mytid)%nlats,nlevp1/)
      allocate(field_reshaped(subdomain_size(1),subdomain_size(2),subdomain_size(3)))
      field_reshaped(:,:,:) = reshape( data_tgcm_intern(:,idx_intern(mytid)%lon0:idx_intern(mytid)%lon1,&
                                                  idx_intern(mytid)%lat0:idx_intern(mytid)%lat1),&
                                subdomain_size, order=(/3,1,2/))
    else if(present(data_state)) then
      subdomain_size = (/idx_nc(mytid)%nlons,idx_nc(mytid)%nlats,size(data_state,dim=1)/)
      allocate(field_reshaped(subdomain_size(1),subdomain_size(2),subdomain_size(3)))
      field_reshaped(:,:,:) = reshape( data_state, &
                                subdomain_size, order=(/3,1,2/))
    else
      write(*,*) "ERROR invalid call to instant_tiegcm_write"
    end if

    ! gather data on root
    if ( mytid == 0) then
      allocate(member(nlon,nlat,size(field_reshaped,dim=3)))
    else
      allocate(member(0,0,0))
    end if

    call gather_model_3d(field_reshaped,&
                         member, &
                         (/idx_nc(mytid)%lon0,idx_nc(mytid)%lat0,levd0/))

    ! write to netcdf
    if ( mytid == 0) then

      start_p(1:3) = 1
      start_p(4) = t_idx

      count_p(1:3)=shape(member)
      count_p(4)=1

      istat = nf90_put_var(ncid, var_id, start=start_p, count=count_p, values=member)
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'instant_tiegcm_write: put '//name)

      istat=nf90_sync(ncid)
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'instant_tiegcm_write: syncronizing file')

    end if

    if (allocated(member)) deallocate(member)
    if (allocated(field_reshaped)) deallocate(field_reshaped)

  end subroutine

  !> Immediately writes one regular-grid field to an open NetCDF file/group, creating dims/vars on demand (for debugging).
  subroutine instant_reg_grid_write(ncid,name,t_idx,data,dataset)
    ! standalone subroutine to write data on tme grid to nc file

    ! tie-gcm
    use mpi_module, only: mytid

    ! intern
    use grid_observation_module, only: reg_grid_dataset_root
    use netcdf_functionality, only: init_temporal_dims, add_model_time
    use quantity_info_module, only: LEVEL_INT,LEVEL_MID

    implicit none

    ! arguments
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer, intent(in) :: t_idx
    real, dimension(:,:,:), intent(in) :: data
    type(reg_grid_dataset_root), intent(in) :: dataset


    ! local
    integer :: istat
    integer :: var_dimensions(4)
    integer :: var_id
    integer :: len_ulim

    integer :: nlon,nlat,nalt

    real, dimension(:,:,:), allocatable :: member

    integer :: start_p(4)
    integer :: count_p(4)

    if ( mytid == 0) then

      ! create dimensions if missing
      istat = nf90_inq_dimid(ncid, "n", var_dimensions(4))
      if (istat /= NF90_NOERR) then
        call init_temporal_dims(ncid)
      end if

      istat = nf90_inq_dimid(ncid, "lon", var_dimensions(1))
      if (istat /= NF90_NOERR) then
        call dataset%add_spatial_dim(ncid)
      end if

      ! get dimension ids
      istat = nf90_inq_dimid(ncid, "lon", var_dimensions(1))
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'instant_reg_grid_write: inquire lon dimension')
      istat = nf90_inq_dimid(ncid, "lat", var_dimensions(2))
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'instant_reg_grid_write: inquire lat dimension')
      istat = nf90_inq_dimid(ncid, "alt", var_dimensions(3))
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'instant_reg_grid_write: inquire alt dimension')

      istat = nf90_inq_dimid(ncid, "n", var_dimensions(4))
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'instant_reg_grid_write: inquire unlimited dimension')

      ! expand time dimension if necessary
      istat = nf90_inquire_dimension(ncid, var_dimensions(4), len=len_ulim)
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'instant_reg_grid_write: inquire length of dimension n')
      if (len_ulim < t_idx) then
        call add_model_time(ncid, t_idx, .false.)
      end if

      ! if variable is not existing create it
      istat = nf90_inq_varid(ncid, trim(name), var_id)
      if (istat /= NF90_NOERR) then
        call rfw_def_var(ncid=ncid,&
                         varname=name,&
                         dimids=var_dimensions,&
                         var_id=var_id)
      end if
    end if

    nlon=size(dataset%lon)
    nlat=size(dataset%lat)
    nalt = (dataset%alt_last-dataset%alt_first)+1

    ! gather data on root
    if ( mytid == 0) then
      allocate(member(nlon,nlat,nalt))
    else
      allocate(member(0,0,0))
    end if

    call gather_model_3d(data,&
                         member, &
                         (/dataset%offset_lon,dataset%offset_lat,1/))

    ! write to netcdf
    if ( mytid == 0) then

      start_p(1:2) = 1
      start_p(3) = dataset%alt_first
      start_p(4) = t_idx

      count_p(1:3)=shape(member)
      count_p(4)=1

      istat = nf90_put_var(ncid, var_id, start=start_p, count=count_p, values=member)
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'instant_reg_grid_write: put '//name)

      istat=nf90_sync(ncid)
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'instant_reg_grid_write: syncronizing file')

    end if

    if (allocated(member)) deallocate(member)

  end subroutine

  !> Computes the expected number of output time steps for the whole run, to preallocate the unlimited dimension.
  function compute_temporal_dimension_size() result (n)

    ! tie-gcm
    use input_module, only: pristart, pristop

    ! intern
    use math_addon_module, only: least_common_multiple
    use configuration, only: cfg_output, cfg_filter

    implicit none

    integer(kind=8),external :: mtime_to_nsec

    integer :: n, lcm

    integer(kind=8) :: nsec_start, nsec_stop

    integer :: delta

    integer :: n_assim_step, n_advance, n_skip

    nsec_start = mtime_to_nsec(pristart(:,1)) ! model start time
    nsec_stop = mtime_to_nsec(pristop(:,1))   ! model stop time

    ! max duration for a model run is a year
    ! a year has 365*86,400 s = 31,536,000 s which is much smaller than max int(kind=4) 2,147,483,647
    ! Thus, it is safe to cast delta to int(kind=4)
    delta = int(nsec_stop-nsec_start,kind=4)

    if(cfg_filter%forecast_duration_sec>0)then
      n_assim_step = (delta-cfg_filter%first_analysis_step_sec)/cfg_filter%forecast_duration_sec + 1
    else
      n_assim_step = 0
    end if

    if(cfg_output%write_every_sec>0)then
      n_advance = delta/cfg_output%write_every_sec + 1
    else
      n_advance = 0
    end if

    if(cfg_filter%forecast_duration_sec>0)then
      lcm = least_common_multiple(cfg_output%write_every_sec,cfg_filter%forecast_duration_sec)
      n_skip = (delta-cfg_filter%first_analysis_step_sec)/lcm + 1
    else
      n_skip = 0
    end if


    n = n_assim_step*2 + n_advance - n_skip

!     write(*,*) "delta ", delta
!     write(*,*) "lcm ", lcm
!     write(*,*) "n_assim_step", n_assim_step
!     write(*,*) "n_advance", n_advance
!     write(*,*) "n_skip", n_skip
!     write(*,*) "n", n

  end function

end module
