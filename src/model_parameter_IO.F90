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
! Writes calibrated model-parameter ensembles/moments to NetCDF and reads perturbation/calibration ensemble files back in.

module model_parameter_IO_module

  type :: model_parameter_writer_type
    integer :: ncid =-1
    integer :: grp_id_members = -1
    integer :: grp_id_moments(4) = -1
    contains
    procedure, pass(this), public :: create => model_parameter_writer_create
    procedure, pass(this), private :: add_variables => model_parameter_writer_add_variables
    procedure, pass(this), public :: sync => model_parameter_writer_sync
    procedure, pass(this), private :: write_param => model_parameter_writer_write_param
    procedure, pass(this), public :: write => model_parameter_writer_write
    procedure, pass(this), public :: close => model_parameter_writer_close

  end type

  type :: model_parameter_reader_type
    integer :: ncid
    integer :: grp_id_perturbations, grp_id_members, grp_id_mean, grp_id_std
    integer :: perturbation_id
    contains
    procedure, pass(this), public :: open => model_parameter_reader_open
    procedure, pass(this), public :: read => model_parameter_reader_read
    procedure, pass(this), public :: close => model_parameter_reader_close
  end type

  integer, parameter, private :: max_dims = 5

  type(model_parameter_writer_type), protected :: model_parameter_writer

  contains


  !> On world rank 0, creates the calibration-parameter NetCDF output file, defines its dimensions/groups, and declares the calibrated parameters' variables.
  subroutine model_parameter_writer_create(this, filename)

    ! extern
    use netcdf

    ! tie-gcm
    USE nchist_module, ONLY: handle_ncerr

    ! intern
    use mod_parallel_pdaf,&
        only: rank_world, n_modeltasks
    use netcdf_functionality,&
        only: add_global_meta_data, init_temporal_dims, init_spacial_dims
    use result_file_writer_module,&
        only: add_member_and_moments_group

    implicit none

    ! arguments
    class(model_parameter_writer_type) :: this
    character(len=*), intent(in) :: filename

    ! local
    integer :: istat

    integer :: dim_id_ulim
    integer :: dim_id_ens
    integer :: dim_id_bin

    if ( rank_world == 0) then

      istat = nf90_create(path=trim(filename),&
                          cmode=NF90_NETCDF4,&
                          ncid=this%ncid)
      if (istat /= NF90_NOERR) call handle_ncerr(istat, 'model_parameter: Error creating '//filename)

      call add_global_meta_data(this%ncid)

      call init_temporal_dims(this%ncid)
      istat = nf90_inquire(this%ncid, unlimiteddimid = dim_id_ulim)

      istat = NF90_DEF_DIM(this%ncid, 'ensemble', n_modeltasks, dim_id_ens)
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error defining dim ensemble')

      istat = NF90_DEF_DIM(this%ncid, 'bin', 37, dim_id_bin)
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error defining dim bin')

      call init_spacial_dims(this%ncid)

      call add_member_and_moments_group(this%ncid,this%grp_id_members,this%grp_id_moments)

      call this%add_variables()
      istat=nf90_sync(this%ncid)

    end if

  end subroutine

  !> On world rank 0, closes this writer's NetCDF file, if open (by model_parameter_writer_create).
  subroutine model_parameter_writer_close(this)

    ! extern
    use netcdf

    ! tie-gcm
    use nchist_module, only: handle_ncerr

    ! intern
    use mod_parallel_pdaf, only: rank_world

    implicit none

    ! arguments
    class(model_parameter_writer_type) :: this

    ! local
    integer :: istat

    if ( rank_world == 0) then
      if (this%ncid>0) then
        istat = nf90_close(this%ncid)
        if (istat /= NF90_NOERR) call handle_ncerr(istat,'model_parameter_writer_close')
      end if
    end if
  end subroutine

  !> Declares a members-group and moments-group NetCDF variable for every calibrated model parameter.
  subroutine model_parameter_writer_add_variables(this)

    ! extern
    use netcdf
    use mpi_f08

    ! intern
    use result_file_writer_module,&
        only: rfw_def_var
    use model_parameter_handling_module,&
        only: model_parameters, PARAM_HANDLING_NONE, PARAM_HANDLING_CALIBRATE

    implicit none

    class(model_parameter_writer_type) :: this

    ! local
    integer, dimension(max_dims) :: var_dimensions
    integer, dimension(max_dims) :: var_dimensions_moments
    integer :: n_nc_dims

    integer :: i,j
    integer :: istat
    integer :: var_id

    do i= 1, size(model_parameters)
      if( model_parameters(i)%handling == PARAM_HANDLING_CALIBRATE ) then

        if(model_parameters(i)%size>1) then
          n_nc_dims=3
          istat = nf90_inq_dimid(this%ncid, trim(model_parameters(i)%nc_dimension), var_dimensions(1))
        else
          n_nc_dims=2
        end if

        istat = nf90_inq_dimid(this%ncid, "ensemble", var_dimensions(n_nc_dims-1))
        istat = nf90_inq_dimid(this%ncid, "n", var_dimensions(n_nc_dims))

        var_dimensions_moments(1:n_nc_dims-2) = var_dimensions(1:n_nc_dims-2)
        var_dimensions_moments(n_nc_dims-1) = var_dimensions(n_nc_dims)


        call rfw_def_var(ncid=this%grp_id_members,&
                         varname=model_parameters(i)%name,&
                         dimids=var_dimensions(1:n_nc_dims),&
                         var_id=var_id)
        istat = nf90_put_att(this%grp_id_members, var_id, &
                             "units", trim(model_parameters(i)%units))
        istat = nf90_put_att(this%grp_id_members, var_id, &
                             "long_name", trim(model_parameters(i)%long_name))

        do j=1,4
          if(this%grp_id_moments(j)>=0)then
            call rfw_def_var(ncid=this%grp_id_moments(j),&
                            varname=model_parameters(i)%name,&
                            dimids=var_dimensions_moments(1:n_nc_dims-1),&
                            var_id=var_id)
            istat = nf90_put_att(this%grp_id_moments(j), var_id, &
                                "units", trim(model_parameters(i)%units))
            istat = nf90_put_att(this%grp_id_moments(j), var_id, &
                                "long_name", trim(model_parameters(i)%long_name))
          end if
        end do

      end if
    end do

  end subroutine

  !> Flushes this writer's NetCDF file to disk.
  subroutine model_parameter_writer_sync(this)

    ! extern
    use netcdf

    implicit none

    class(model_parameter_writer_type) :: this
    integer :: istat

    istat=nf90_sync(this%ncid)

  end subroutine

  !> If the handling of cal_par is set up to be calibrated: collects this parameter's
  !! ensemble of calibrated values (distributed in state_cal over all model instances)
  !! on world rank 0, computes the moments, and writes them together with the ensemble
  !! to the NetCDF file.
  subroutine model_parameter_writer_write_param(this,cal_par,state_cal)

    ! extern
    use mpi_f08
    use netcdf

    ! tie-gcm
    use mpi_module, only: handle_mpi_err
    use nchist_module, only: handle_ncerr

    ! intern
    use array_print_module,&
        only: printMat
    use array_mapping_module,&
        only: flatten
    use model_parameter_handling_module,&
        only: model_parameter, PARAM_HANDLING_NONE, PARAM_HANDLING_CALIBRATE
    use mod_assimilation,&
        only: analysis_step_dynamics_count
    use mod_parallel_pdaf,&
        only: rank_world, n_modeltasks, COMM_Couple
    use netcdf_functionality,&
        only: add_model_time
    use PDAF,&
        only: PDAF_diag_compute_moments

    implicit none

    ! arguments
    class(model_parameter_writer_type) :: this
    type(model_parameter) :: cal_par
    real, dimension(:), intent(in) :: state_cal

    ! local
    integer, dimension(max_dims) :: start_p
    integer, dimension(max_dims) :: count_p

    integer :: istat, ierr
    integer :: i,j
    integer :: var_id

    integer :: ndims
    integer, dimension(max_dims) :: dimids

    real, dimension(:,:), allocatable :: ensemble, moments

    if( cal_par%handling == PARAM_HANDLING_CALIBRATE ) then

        if(rank_world==0)then
          allocate(ensemble(size(state_cal),n_modeltasks))
        else
          allocate(ensemble(0,0))
        end if

        call MPI_Gather(sendbuf=state_cal,&
                sendcount=size(state_cal),&
                sendtype=MPI_REAL8,&
                recvbuf=ensemble,&
                recvcount=size(state_cal),&
                recvtype=MPI_REAL8,&
                root=0,&
                comm=COMM_Couple,&
                ierror=ierr)
        if (ierr /= 0) call handle_mpi_err(ierr,'model_parameter_write MPI_Gather')

        if(rank_world==0)then

          allocate(moments(size(state_cal),4))

          call PDAF_diag_compute_moments( dim_p=size(state_cal),&
                                          dim_ens=n_modeltasks,&
                                          ens=ensemble,&
                                          kmax=4,&
                                          moments=moments)
          moments(:,2) = sqrt(moments(:,2))

          istat = nf90_inq_varid(this%grp_id_members,cal_par%name,var_id)
          if (istat /= NF90_NOERR) call handle_ncerr(istat,'model_parameter_writer_write')
          istat = nf90_inquire_variable(ncid=this%grp_id_members, &
                                        varid=var_id, &
                                        ndims=ndims, &
                                        dimids=dimids)
          if (istat /= NF90_NOERR) call handle_ncerr(istat,'model_parameter_writer_write')

          istat=nf90_inquire_dimension(this%grp_id_members,dimids(ndims),len=start_p(ndims) )
          if (istat /= NF90_NOERR) call handle_ncerr(istat,'model_parameter_writer_write')

          ! first entry is inital ensemble: Add one to analysis_step_count
          if (start_p(ndims) < analysis_step_dynamics_count+1) then
            call add_model_time(this%ncid, analysis_step_dynamics_count+1, .false.)
          end if

          start_p(1:ndims-1) = 1
          start_p(ndims) = analysis_step_dynamics_count+1 ! starts at zero


          do i=1,ndims-1
            istat=nf90_inquire_dimension(this%grp_id_members,dimids(i),len=count_p(i))
             if (istat /= NF90_NOERR) call handle_ncerr(istat,'model_parameter_writer_write')
          end do
          count_p(ndims)=1

          istat = nf90_put_var(this%grp_id_members, &
                                var_id, &
                                start=start_p(1:ndims), &
                                count=count_p(1:ndims), &
                                values=ensemble)
          if (istat /= NF90_NOERR) call handle_ncerr(istat,'model_parameter_writer_write: put '//cal_par%name)

          ! remove entry for ensemble
          start_p(ndims-1) = start_p(ndims)
          count_p(ndims-1) = count_p(ndims)
          ndims = ndims-1

          do j=1,4
            if(this%grp_id_moments(j)>=0)then
              istat = nf90_inq_varid(this%grp_id_moments(j),cal_par%name,var_id)
              istat = nf90_put_var(this%grp_id_moments(j), &
                            var_id, &
                            start=start_p(1:ndims), &
                            count=count_p(1:ndims), &
                            values=moments(:,j))
            end if
          end do

          deallocate(moments)

        end if
        deallocate(ensemble)

    end if


  end subroutine

  !> Loops over all calibrated parameters and writes them to the netcdf file calling model_parameter_writer_write_param.
  subroutine model_parameter_writer_write(this,state_p)

    ! tie-gcm
    use mpi_module, only: mytid

    ! intern
    use model_parameter_handling_module, only: model_parameters
    use state_module, only: state_vector

    implicit none

    class(model_parameter_writer_type) :: this
    REAL, INTENT(inout) :: state_p(:)  ! local state vector

    integer :: i
    integer :: idx

    if(mytid==0)then
      do i= state_vector%idx_cal_0, state_vector%idx_cal_1

          idx = findloc( model_parameters%name,state_vector%map%fd_name(i),dim=1 )

          call this%write_param( &
                  model_parameters(idx), &
                  state_p(state_vector%map%idx_R(i, 0)%begin_p : &
                          state_vector%map%idx_R(i, 0)%back_p))

      end do
      call this%sync()
    end if

  end subroutine


  !> Opens an ensemble/perturbation NetCDF file, locates its known groups, validates it has enough members, and determines which ensemble member this task should read (honoring an optional skip_list).
  subroutine model_parameter_reader_open(this, ensemble_file, skip_list)

    ! extern
    use mpi_f08
    use netcdf

    ! tie-gcm
    use nchist_module, only: handle_ncerr

    ! intern
    use mod_parallel_pdaf, only: n_modeltasks, task_id
    use result_file_writer_module, only: set_collective_access_of_all_vars


    implicit none

    ! arguments
    class(model_parameter_reader_type) :: this
    character(len=*), intent(in) :: ensemble_file
    integer, dimension(:), intent(in), optional :: skip_list

    ! local
    integer :: istat
    integer :: dim_id_ens, dim_len_ens
    character(len=256) :: dskfile

    integer :: pos
    integer :: cumsum
    logical, dimension(:), allocatable :: mask
    logical :: use_skip_list
    character(len=256) :: errmsg

    call getfile(ensemble_file,dskfile)

    istat = nf90_open(path=trim(dskfile), &
                      mode=NF90_NOWRITE, &
                      ncid=this%ncid)!, &
                      !comm=MPI_COMM_WORLD%MPI_VAL, &
                      !info=MPI_INFO_NULL%MPI_VAL)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'read_calibration_ensemble')

    istat = nf90_inq_ncid(this%ncid, 'perturbations', this%grp_id_perturbations)
    if (istat /= NF90_NOERR) write(*,*) 'WARNING netcdf file does not contain group perturbations'

    istat = nf90_inq_ncid(this%ncid, 'members', this%grp_id_members)
    if (istat /= NF90_NOERR) write(*,*) 'WARNING netcdf file does not contain group members'

    istat = nf90_inq_ncid(this%ncid, 'mean', this%grp_id_mean)
    if (istat /= NF90_NOERR) write(*,*) 'WARNING netcdf file does not contain group mean'

    istat = nf90_inq_ncid(this%ncid, 'std', this%grp_id_std)
    if (istat /= NF90_NOERR) write(*,*) 'WARNING netcdf file does not contain group std'

    ! independent access is extremly slow on cluster noctua2
    ! call set_collective_access_of_all_vars(this%grp_id_perturbations)
    ! call set_collective_access_of_all_vars(this%grp_id_members)

    dim_id_ens = -1
    istat = nf90_inq_dimid(this%ncid, "ensemble", dim_id_ens)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'read_calibration_ensemble: finding dimension ensemble')
    dim_len_ens = 0
    istat = nf90_inquire_dimension(ncid=this%ncid, dimid=dim_id_ens, len=dim_len_ens)
    if(dim_len_ens<n_modeltasks) then
      call shutdown('number of ensemble members in '//trim(ensemble_file)//' is samller than members required for this run')
    end if

    ! define perturbation_id used for this task
    use_skip_list = .false.
    if (present(skip_list)) use_skip_list = size(skip_list) > 0

    if (use_skip_list) then

      if(n_modeltasks+size(skip_list)>dim_len_ens)then
        write(errmsg,*) "The perturbation file does not contain enough members for the ensemble!",&
        "the perturbation files contains ", dim_len_ens, " members, ", &
        size(skip_list), " perturbation members are skipped. ", &
        n_modeltasks, " are required"
        call shutdown(trim(errmsg))
      end if

      allocate(mask(n_modeltasks+size(skip_list)))
      mask = .true.
      mask(skip_list) = .false.
      cumsum=0
      do pos=0, size(mask)
        if (mask(pos)) cumsum=cumsum+1
        if(cumsum==task_id) exit
      end do
      this%perturbation_id = pos
      write(*,*) "reading perturbation member",  this%perturbation_id, &
        " for ensemble member ", task_id

      deallocate(mask)
    else
      this%perturbation_id = task_id
    end if

  end subroutine

  !> Closes this reader's NetCDF file.
  subroutine model_parameter_reader_close(this)

    ! extern
    use netcdf

    ! tie-gcm
    USE nchist_module, ONLY: handle_ncerr

    implicit none

    ! arguments
    class(model_parameter_reader_type) :: this

    ! local
    integer :: istat

    istat = nf90_close(this%ncid)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'model_parameter_reader_close')

  end subroutine

  !> If m_parameter has a handling mode, reads this task's ensemble member's value(s) from the appropriate NetCDF group, either as a constant or as a time-interpolator built from the time-variable data.
  subroutine model_parameter_reader_read(this, m_parameter)

    ! intern
    use netcdf_functionality, only: read_time_variable_and_convert
    use model_parameter_handling_module,&
      only: model_parameter,&
            PARAM_HANDLING_NONE,&
            PARAM_HANDLING_PERTURB,&
            PARAM_HANDLING_CALIBRATE,&
            PARAM_HANDLING_OVERWRITE,&
            PARAM_HANDLING_OVERWRITE_MEAN

    ! apmg_fortran_lib
    use search_module, only: NN_LEFT

    implicit none

    ! arguments
    class(model_parameter_reader_type) :: this
    type(model_parameter), intent(inout) :: m_parameter

    ! local
    real, dimension(:), allocatable :: time
    real, dimension(:,:), allocatable :: values

    character(len=16) :: time_dim_name

    integer :: j

    if(m_parameter%handling/=PARAM_HANDLING_NONE) then

      select case(m_parameter%handling)
        case(PARAM_HANDLING_PERTURB)
          call get_time_dim_name(this%grp_id_perturbations,&
                                 m_parameter%name,&
                                 time_dim_name,&
                                 m_parameter%time_variable_perturbation)
        case(PARAM_HANDLING_CALIBRATE,PARAM_HANDLING_OVERWRITE)
          call get_time_dim_name(this%grp_id_members,&
                                 m_parameter%name,&
                                 time_dim_name,&
                                 m_parameter%time_variable_perturbation)
        case(PARAM_HANDLING_OVERWRITE_MEAN)
          call get_time_dim_name(this%grp_id_mean,&
                                 m_parameter%name,&
                                 time_dim_name,&
                                 m_parameter%time_variable_perturbation)
      end select

      if(m_parameter%time_variable_perturbation)then
        write(*,*) '"', trim(m_parameter%name), '" has time variable perturbations'
      else
        write(*,*) '"', trim(m_parameter%name), '" has constant perturbation'
      end if

      ! time in seconds after model start
      if(m_parameter%time_variable_perturbation) then
        call read_time_variable_and_convert(this%ncid,time,time_dim_name)
      end if

      if(m_parameter%time_variable_perturbation) then
        allocate(values(m_parameter%size,size(time)))
      else
        allocate(values(m_parameter%size,1))
      end if

      select case(m_parameter%handling)
        case(PARAM_HANDLING_PERTURB)
          call read_var_in_group(this%grp_id_perturbations, m_parameter%name, values, this%perturbation_id )
        case(PARAM_HANDLING_CALIBRATE,PARAM_HANDLING_OVERWRITE)
          call read_var_in_group(this%grp_id_members, m_parameter%name, values, this%perturbation_id )
        case(PARAM_HANDLING_OVERWRITE_MEAN)
          call read_var_in_group(this%grp_id_mean, m_parameter%name, values, this%perturbation_id )
      end select

      if(m_parameter%time_variable_perturbation) then
        do j=1,m_parameter%size
          call m_parameter%interp(j)%init(x_src=time,&
                                          y_src=values(j,:),&
                                          degree=3)
                                          !neighbour=NN_LEFT,&
                                          !assume_sorted=.true.)
        end do
      else
        m_parameter%val = values(:,1)
      end if

      deallocate(values)
      if(allocated(time)) deallocate(time)

    end if

  end subroutine

  !> Finds the given variable's time-like dimension (name containing "time", falling back to "n") and reports whether it has more than one entry (i.e. is time-variable).
  subroutine get_time_dim_name(grp_id,name,time_dim_name,is_timevariable)

    ! extern
    use netcdf

    ! tie-gcm
    use nchist_module, only: handle_ncerr

    ! intern
    use character_routines_module, only: string_contains

    implicit none

    ! arguments
    integer, intent(in) :: grp_id
    character(len=*), intent(in) :: name

    character(len=16), intent(out) :: time_dim_name
    logical, intent(out) :: is_timevariable

    ! local
    integer :: istat
    integer :: var_id

    integer :: ndims
    integer :: dimids(5) ! we do not expect higher dimensional variables

    integer :: dim_time_len, dim_len
    character(len=16) :: dim_name

    integer :: i

    istat = nf90_inq_varid(ncid=grp_id, name=name, varid=var_id)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'inquire id '//name)

    istat = nf90_inquire_variable(ncid=grp_id, &
                                  varid=var_id, &
                                  ndims=ndims)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'get_time_dim_name')

    istat = nf90_inquire_variable(ncid=grp_id, &
                                  varid=var_id, &
                                  dimids=dimids)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'get_time_dim_name')

    time_dim_name = 'time' ! default value
    dim_time_len = 0
    do i=1,ndims
      istat = nf90_inquire_dimension(ncid=grp_id, dimid=dimids(i), name=dim_name, len=dim_len)
      if(string_contains(dim_name,'time'))then
        dim_time_len = dim_len
        time_dim_name = dim_name
        exit
      else if(dim_name=='n')then
        write(*,*) 'WARNING found dimension n instead of time. Looking further for time.'
        dim_time_len = dim_len
      end if
    end do


    if(dim_time_len>1)then
      is_timevariable = .true.
    else
      is_timevariable = .false.
    end if

  end subroutine

  !> Reads the named NetCDF variable into values. The start position in values is determined by perturbation_id.
  subroutine read_var_in_group(grp_id,name,values,perturbation_id)

    ! extern
    use netcdf

    ! tie-gcm
    use nchist_module, only: handle_ncerr

    ! intern

    ! apmg_fortran_lib
    use array_print_module, only: printMat

    implicit none

    ! arguments
    integer, intent(in) :: grp_id
    character(len=*), intent(in) :: name
    real, dimension(:,:), intent(out) :: values
    integer, intent(in) :: perturbation_id

    ! local
    integer :: istat
    integer :: var_id

    integer :: ndims

    integer, dimension(:), allocatable :: dimids
    integer, dimension(:), allocatable :: start_p
    integer, dimension(:), allocatable :: count_p

    character(len=32) :: dim_name
    integer :: dim_len

    integer :: i

    istat = nf90_inq_varid(ncid=grp_id, name=name, varid=var_id)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'inquire id '//name)


    istat = nf90_inquire_variable(ncid=grp_id, &
                                  varid=var_id, &
                                  ndims=ndims)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'model_parameter_reader_read')

    allocate(dimids(ndims))

    istat = nf90_inquire_variable(ncid=grp_id, &
                                  varid=var_id, &
                                  dimids=dimids)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'model_parameter_reader_read')

    allocate(start_p(ndims))
    allocate(count_p(ndims))

    start_p = 1

    ! ndims   : time
    ! ndims-1 : ensemble
    ! ndims-2 : other

    do i=1,ndims
      istat=nf90_inquire_dimension(grp_id,dimids(i),len=dim_len,name=dim_name )
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'model_parameter_reader_read')
      select case(trim(dim_name))
        case('ensemble')
          count_p(i)=1
          start_p(i)=perturbation_id
        case default
          count_p(i)=dim_len
      end select

    end do

!       write(*,*) 'count_p: ', count_p
!       write(*,*) 'start_p: ', start_p

    istat = nf90_get_var(ncid=grp_id, &
                         varid=var_id, &
                         values=values, &
                         start=start_p, &
                         count=count_p)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'model_parameter_reader_read')
!       call printMat(values)

    deallocate(dimids)
    deallocate(start_p)
    deallocate(count_p)

  end subroutine


end module model_parameter_IO_module
