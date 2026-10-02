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
! Shared NetCDF helper routines: spatial/temporal dimension setup, time read/write, group handling, and metadata for output files.

module netcdf_functionality

  ! extern
  use netcdf

  ! tie-gcm
  use nchist_module, only: handle_ncerr

  implicit none

  character(len=33), protected :: start_time

  contains

  !> Broadcasts the run's creation date/time (determined on world rank 0) to all MPI ranks for use in output metadata.
  subroutine init_start_time

    ! extern
    use mpi_f08

    ! intern
    use mod_parallel_pdaf, only: rank_world

    implicit none
    character(len=16) :: create_date,create_time
    integer :: ierr

    if (rank_world==0) then
      ! datetime is defined in util.F
      call datetime(create_date,create_time)
      start_time = trim(create_date)//' '//trim(create_time)
    end if

    call MPI_Bcast( start_time, 33, MPI_CHARACTER, 0, &
               MPI_COMM_WORLD, ierr )

  end subroutine

  !> Defines the spatial dimensions and coordinate variables for the sundomain of
  !! the calling rank only, plus a time level dimension.
  !! dim_ids and var_ids are returned in the order (lon, lat, lev, ilev, timelevel).
  subroutine init_spacial_dims_subdomain(ncid, nc_lon0, nc_lon1, nc_lat0, nc_lat1, &
                                         dim_ids, var_ids)

    ! tie-gcm
    use nchist_module, only: handle_ncerr
    use params_module, only: nlevp1

    implicit none

    ! arguments
    integer, intent(in) :: ncid
    integer, intent(in) :: nc_lon0, nc_lon1, nc_lat0, nc_lat1
    integer, intent(out) :: dim_ids(5)
    integer, intent(out) :: var_ids(5)

    ! local
    integer :: istat

    ! Define dimensions --------------------------------------

    istat = nf90_def_dim(ncid, 'lon', nc_lon1-nc_lon0+1, dim_ids(1))
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error defining dim lon')

    istat = nf90_def_dim(ncid, 'lat', nc_lat1-nc_lat0+1, dim_ids(2))
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error defining dim lat')

    istat = nf90_def_dim(ncid, 'lev', nlevp1, dim_ids(3))
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error defining dim lev')

    istat = nf90_def_dim(ncid, 'ilev', nlevp1, dim_ids(4))
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error defining dim ilev')

    ! itc and itp swap every model step (see advance.F) The caller therefore writes
    ! itp into index 1 and itc into index 2
    istat = nf90_def_dim(ncid, 'timelevel', 2, dim_ids(5))
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error defining dim timelevel')

    ! define variables -------------------------------------------------

    ! longitude
    istat = nf90_def_var(ncid,'lon', NF90_DOUBLE, dim_ids(1), var_ids(1))
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error defining longitude dimension variable')

    istat = nf90_put_att(ncid,var_ids(1),"long_name", "geographic longitude (-west, +east)")
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error putting attribute to lon')
    istat = nf90_put_att(ncid,var_ids(1),"units",'degrees_east')
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error putting attribute to lon')

    ! latitude
    istat = nf90_def_var(ncid,'lat', NF90_DOUBLE, dim_ids(2), var_ids(2))
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error defining latitude dimension variable')

    istat = nf90_put_att(ncid,var_ids(2),"long_name", "geographic latitude (-south, +north)")
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error putting attribute to lat')
    istat = nf90_put_att(ncid,var_ids(2),"units",'degrees_north')
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error putting attribute to lat')

    ! Midpoint levels coordinate variable lev(lev):
    istat = nf90_def_var(ncid,"lev",NF90_DOUBLE, dim_ids(3), var_ids(3))
    if (istat /= NF90_NOERR) call handle_ncerr(istat, 'Error defining midpoint levels coordinate variable (lev)')
    istat = nf90_put_att(ncid,var_ids(3),"long_name", "midpoint levels")
    istat = nf90_put_att(ncid,var_ids(3),"short_name", "ln(p0/p)")
    istat = nf90_put_att(ncid,var_ids(3),"units"," ")
    istat = nf90_put_att(ncid,var_ids(3),"positive",'up')
    istat = nf90_put_att(ncid,var_ids(3),"standard_name", "atmosphere_ln_pressure_coordinate")
    istat = nf90_put_att(ncid,var_ids(3),"formula_terms", "p0: p0 lev: lev")
    istat = nf90_put_att(ncid,var_ids(3),"formula", "p(k) = p0 * exp(-lev(k))")

    ! Interface levels coordinate array:
    istat = nf90_def_var(ncid,"ilev",NF90_DOUBLE, dim_ids(4), var_ids(4))
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error defining interface levels coordinate variable ilev')
    istat = nf90_put_att(ncid,var_ids(4),"long_name", "interface levels")
    istat = nf90_put_att(ncid,var_ids(4),"short_name", "ln(p0/p)")
    istat = nf90_put_att(ncid,var_ids(4),"units"," ")
    istat = nf90_put_att(ncid,var_ids(4),"positive", "up")
    istat = nf90_put_att(ncid,var_ids(4),"standard_name", "atmosphere_ln_pressure_coordinate")
    istat = nf90_put_att(ncid,var_ids(4),"formula_terms", "p0: p0 lev: ilev")
    istat = nf90_put_att(ncid,var_ids(4),"formula", "p(k) = p0 * exp(-ilev(k))")

    ! time level
    istat = nf90_def_var(ncid,"timelevel",NF90_INT, dim_ids(5), var_ids(5))
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error defining timelevel coordinate variable')
    istat = nf90_put_att(ncid,var_ids(5),"long_name", "leapfrog time level")
    istat = nf90_put_att(ncid,var_ids(5),"units"," ")
    istat = nf90_put_att(ncid,var_ids(5),"description", &
      "1: previous time level (itp), 2: current time level (itc)")

  end subroutine init_spacial_dims_subdomain

  !> Writes the coordinate values belonging to init_spacial_dims_subdomain.
  !! Must be called AFTER the caller has left define mode (NF90_ENDDEF).
  subroutine write_spacial_dims_subdomain(ncid, nc_lon0, nc_lon1, nc_lat0, nc_lat1, var_ids)

    ! tie-gcm
    use nchist_module, only: handle_ncerr
    use params_module, only: glon, glat, zpmid, zpint

    implicit none

    ! arguments
    integer, intent(in) :: ncid
    integer, intent(in) :: nc_lon0, nc_lon1, nc_lat0, nc_lat1
    integer, intent(in) :: var_ids(5)

    ! local
    integer :: istat

    istat = nf90_put_var(ncid, var_ids(1), glon(nc_lon0:nc_lon1))
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error giving values to glon coord var')
    istat = nf90_put_var(ncid, var_ids(2), glat(nc_lat0:nc_lat1))
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error giving values to glat coord var')

    istat = nf90_put_var(ncid, var_ids(3), zpmid)  ! midpoint levels
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error giving values to zpmid coord var')
    istat = nf90_put_var(ncid, var_ids(4), zpint) ! interface levels
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error giving values to zpint coord var')

    istat = nf90_put_var(ncid, var_ids(5), (/1,2/))
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error giving values to timelevel coord var')

  end subroutine write_spacial_dims_subdomain

  !> Defines the lon/lat/lev/ilev spatial dimensions and coordinate variables in a NetCDF file and writes their values.
  subroutine init_spacial_dims(ncid)

    ! tie-gcm
    use nchist_module, only: handle_ncerr
    use params_module, only: nlon, nlat, glon, glat, zpmid, zpint, nlevp1

    implicit none

    ! arguments
    integer, intent(in) :: ncid

    ! local
    integer :: istat
    integer :: var_id_lat, var_id_lon, var_id_ilev, var_id_lev ! coordinate variables
    integer :: dim_id_lon, dim_id_lat,  dim_id_ilev, dim_id_lev, dim_id_lev1

    ! Define dimensions --------------------------------------

    istat = nf90_def_dim(ncid, 'lon', nlon, dim_id_lon)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error defining dim lon')

    istat = nf90_def_dim(ncid, 'lat', nlat, dim_id_lat)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error defining dim lat')

    istat = nf90_def_dim(ncid, 'ilev', nlevp1, dim_id_ilev)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error defining dim ilev')

    istat = nf90_def_dim(ncid, 'lev', nlevp1, dim_id_lev)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error defining dim lev')

    ! Quantities without vertical extent (LEVEL_NONE), e.g. VTEC, are written
    ! with this degenerate vertical dimension. This keeps the rank of all
    ! variables identical, so that the gathering and writing routines do not
    ! have to distinguish between two and three dimensional quantities.
    istat = nf90_def_dim(ncid, 'lev1', 1, dim_id_lev1)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error defining dim lev1')

    ! define variables -------------------------------------------------

    ! longitude
    istat = nf90_def_var(ncid,'lon', NF90_DOUBLE, dim_id_lon, var_id_lon)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error defining longitude dimension variable')


    istat = nf90_put_att(ncid,var_id_lon,"long_name", "geographic longitude (-west, +east)")
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error putting attribute to lon')
    istat = nf90_put_att(ncid,var_id_lon,"units",'degrees_east')
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error putting attribute to lon')

    ! latitude
    istat = nf90_def_var(ncid,'lat', NF90_DOUBLE, dim_id_lat, var_id_lat)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error defining latitude dimension variable')

    istat = nf90_put_att(ncid,var_id_lat,"long_name", "geographic latitude (-south, +north)")
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error putting attribute to lat')
    istat = nf90_put_att(ncid,var_id_lat,"units",'degrees_north')
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error putting attribute to lat')

    ! Interface levels coordinate array:
    istat = nf90_def_var(ncid,"ilev",NF90_DOUBLE, dim_id_ilev, var_id_ilev)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error defining interface levels coordinate variable ilev')
    istat = nf90_put_att(ncid,var_id_ilev,"long_name", "interface levels")
    istat = nf90_put_att(ncid,var_id_ilev,"short_name", "ln(p0/p)")
    istat = nf90_put_att(ncid,var_id_ilev,"units"," ")
    istat = nf90_put_att(ncid,var_id_ilev,"positive", "up")
    istat = nf90_put_att(ncid,var_id_ilev,"standard_name", "atmosphere_ln_pressure_coordinate")
    istat = nf90_put_att(ncid,var_id_ilev,"formula_terms", "p0: p0 lev: ilev")
    istat = nf90_put_att(ncid,var_id_ilev,"formula", "p(k) = p0 * exp(-ilev(k))")

    ! Midpoint levels coordinate variable lev(lev):
    istat = nf90_def_var(ncid,"lev",NF90_DOUBLE, dim_id_lev, var_id_lev)
    if (istat /= NF90_NOERR) call handle_ncerr(istat, 'Error defining midpoint levels coordinate variable (lev)')
    ! long name of lev coord var:
    istat = nf90_put_att(ncid,var_id_lev,"long_name", "midpoint levels")
    istat = nf90_put_att(ncid,var_id_lev,"short_name", "ln(p0/p)")
    istat = nf90_put_att(ncid,var_id_lev,"units"," ")
    istat = nf90_put_att(ncid,var_id_lev,"positive",'up')
    istat = nf90_put_att(ncid,var_id_lev,"standard_name", "atmosphere_ln_pressure_coordinate")
    istat = nf90_put_att(ncid,var_id_lev,"formula_terms", "p0: p0 lev: lev")
    istat = nf90_put_att(ncid,var_id_lev,"formula", "p(k) = p0 * exp(-lev(k))")

    ! end of definition --------------------------------------------

    istat = NF90_ENDDEF(ncid)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error leaving define mode')

    istat = nf90_put_var(ncid,var_id_lon,glon)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error giving values to glon coord var')
    istat = nf90_put_var(ncid,var_id_lat,glat)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error giving values to glat coord var')

    istat = nf90_put_var(ncid, var_id_lev, zpmid )  ! midpoint levels
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error giving values to zpmid coord var')
    istat = nf90_put_var(ncid, var_id_ilev, zpint ) ! interface levels
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error giving values to zpint coord var')
  end subroutine init_spacial_dims

  !> Defines the time-related dimensions and variables (mtime, time, doy, model_step) in a NetCDF file.
  subroutine init_temporal_dims(ncid, dim_t)

    ! tie-gcm
    use hist_module, only: h ,sh
    use nchist_module, only: handle_ncerr

    implicit none

    ! arguments
    integer, intent(in) :: ncid
    integer, intent(in), optional :: dim_t

    ! local
    integer :: istat
    character(len=80) :: char80

    integer :: imo,ida,startmtime(4)

    integer :: dim_id_mtime, dim_id_ulim

    integer :: var_id_doy
    integer :: var_id_time
    integer :: var_id_mtime
    integer :: var_id_step

    real :: rmins
    real,external :: mtime_to_datestr

      istat = nf90_def_dim(ncid, 'mtimedim', 3, dim_id_mtime)
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error defining dim mtime')

      if(present(dim_t))then
        istat = nf90_def_dim(ncid, "n", dim_t, dim_id_ulim)
        if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error defining n dimension')
      else
        istat = nf90_def_dim(ncid, "n", NF90_UNLIMITED, dim_id_ulim)
        if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error defining n dimension')
      end if

      ! define time variables ------------------------------------------------

      ! mtime
      istat = nf90_def_var(ncid,"mtime",NF90_INT, (/dim_id_mtime, dim_id_ulim/), var_id_mtime)
      if (istat /= NF90_NOERR) call handle_ncerr(istat, 'Error defining mtime variable')
      istat = nf90_put_att(ncid, var_id_mtime, "long_name", "model times (doy, hour, minute)")
      if (istat /= NF90_NOERR) call handle_ncerr(istat, 'Error defining long_name attribute of mtime variable')
      istat = nf90_put_att(ncid, var_id_mtime,"units", "day, hour, minute")
      if (istat /= NF90_NOERR) call handle_ncerr(istat, 'Error defining units of mtime variable')

      ! Time (coordinate variable time(time)). This is days since
      ! the initial run's start time. The units string is: yyyy-m-d,
      ! where yyyy is the year, m is month, and d is day of the source
      ! start time.
      !
      istat = nf90_def_var(ncid,"time",NF90_DOUBLE, dim_id_ulim, var_id_time)
      if (istat /= NF90_NOERR) call handle_ncerr(istat, 'Error defining time dimension variable')

      ! day of year
      istat = nf90_def_var(ncid,"doy",NF90_DOUBLE, dim_id_ulim, var_id_doy)
      if (istat /= NF90_NOERR) call handle_ncerr(istat, 'Error defining doy variable')
      istat = nf90_put_att(ncid, var_id_doy,"long_name", "day of year including fraction of day")
      istat = nf90_put_att(ncid, var_id_mtime,"units", "days")


      startmtime(1:3) = sh%initial_mtime(1:3) ; startmtime(4) = 0

      rmins = mtime_to_datestr(sh%initial_year,startmtime,imo,ida,char80)

      istat = nf90_put_att(ncid, var_id_time,"long_name","time")
      istat = nf90_put_att(ncid, var_id_time,"units", trim(char80))
      istat = nf90_put_att(ncid, var_id_time,"initial_year", h%initial_year)
      istat = nf90_put_att(ncid, var_id_time,"initial_day", h%initial_day)
      istat = nf90_put_att(ncid, var_id_time,"initial_mtime", h%initial_mtime)

      ! model step
      istat = nf90_def_var(ncid,"model_step",NF90_INT, dim_id_ulim, var_id_step)
      if (istat /= NF90_NOERR) call handle_ncerr(istat, 'Error defining model step dimension variable')


      istat = nf90_put_att(ncid, var_id_step, "long_name", "model steps after advancing inital values")
      istat = nf90_put_att(ncid, var_id_step,"units", "none")

  end subroutine init_temporal_dims

  !> Writes the current model time (mtime, time, model_step, doy) to record counter of a NetCDF file.
  subroutine add_model_time(ncid, counter, collective)

    ! tie-gcm
    use hist_module, only: modeltime, sh
    use init_module, only: istep

    implicit none

    ! arguments
    integer, intent(in) :: ncid
    integer, intent(in) :: counter
    logical, intent(in) :: collective

    ! local
    integer :: istat
    integer :: mtimeinit(4)

    real :: rmins
    real,external :: mtime_delta

    integer :: var_id_doy
    integer :: var_id_time
    integer :: var_id_mtime
    integer :: var_id_step

    istat = nf90_inq_varid(ncid, "mtime", var_id_mtime)
    if (istat /= NF90_NOERR) call handle_ncerr(istat, 'Error INQUIRE mtime')

    istat = nf90_inq_varid(ncid, "time", var_id_time)
    if (istat /= NF90_NOERR) call handle_ncerr(istat, 'Error INQUIRE time')

    istat = nf90_inq_varid(ncid, "model_step", var_id_step)
    if (istat /= NF90_NOERR) call handle_ncerr(istat, 'Error INQUIRE model_step')

    istat = nf90_inq_varid(ncid, "doy", var_id_doy)
    if (istat /= NF90_NOERR) call handle_ncerr(istat, 'Error INQUIRE doy')

    if(collective) then
      istat = nf90_var_par_access(ncid, var_id_mtime, NF90_COLLECTIVE);
      if (istat /= NF90_NOERR) call handle_ncerr(istat, 'Error setting par_access mtime')

      istat = nf90_var_par_access(ncid, var_id_time, NF90_COLLECTIVE);
      if (istat /= NF90_NOERR) call handle_ncerr(istat, 'Error setting par_access time')

      istat = nf90_var_par_access(ncid, var_id_step, NF90_COLLECTIVE);
      if (istat /= NF90_NOERR) call handle_ncerr(istat, 'Error setting par_access model_step')

      istat = nf90_var_par_access(ncid, var_id_doy, NF90_COLLECTIVE);
      if (istat /= NF90_NOERR) call handle_ncerr(istat, 'Error setting par_access doy')
    end if

    istat = nf90_put_var(ncid, var_id_mtime, start=(/1,counter/), count=(/3, 1/),  values=modeltime(1:3))
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error putting mtime')

    mtimeinit(1:3) = sh%initial_mtime(1:3) ; mtimeinit(4) = 0
    rmins = mtime_delta(mtimeinit, modeltime)

    istat = nf90_put_var(ncid=ncid, varid=var_id_time, values=rmins, start=(/counter/))
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error putting time coordinate')

    istat = nf90_put_var(ncid, var_id_step, values=istep, start=(/counter/))
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error putting model_step')


    istat = nf90_put_var(ncid, var_id_doy, start=(/counter/), &
      values=modeltime(1) + modeltime(2)/24. + modeltime(3)/1440. +  modeltime(4)/86400.)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error putting doy')

  end subroutine

  !> Writes the current wall-clock date/time as a global "create_date" attribute.
  subroutine add_creation_time(ncid)

    implicit none
    ! arguments
    integer, intent(in) :: ncid

    ! local
    integer :: istat
    character(len=80) :: char80
    character(len=24) :: create_date,create_time

    ! datetime is defined in util.F
    call datetime(create_date,create_time)
    char80 = trim(create_date)//' '//trim(create_time)
    istat = nf90_put_att(ncid,NF90_GLOBAL, "create_date", trim(char80))
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error writing create time')

  end subroutine

  !> Writes global metadata attributes (run start time, creation time, title, repository URL, git versions of tiegcm-pdaf/pdaf-binding-tiegcm/pdaf/geodetic-fortran-utilities/esmf) to a NetCDF file.
  subroutine add_global_meta_data(ncid)

    implicit none
#include "gitversion.inc"

    ! arguments
    integer, intent(in) :: ncid

    integer :: istat

    istat = nf90_put_att(ncid, NF90_GLOBAL, 'run_start_time', start_time)
    call add_creation_time(ncid)

    istat = nf90_put_att(ncid, NF90_GLOBAL, 'title', 'TIE-GCM 3.0 PDAF output')
    istat = nf90_put_att(ncid, NF90_GLOBAL, 'software repository', 'https://github.com/rainbowsend/tiegcm')

    ! gitversion_* variables are in gitversion.inc, created by the makefile
    ! from the checked-out commit of the tiegcm-pdaf repo and each deps/
    ! submodule.
    istat = nf90_put_att(ncid, NF90_GLOBAL, 'git version tiegcm-pdaf', gitversion_tiegcm_pdaf)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error writing git version tiegcm-pdaf attribute')
    istat = nf90_put_att(ncid, NF90_GLOBAL, 'git version pdaf-binding-tiegcm', gitversion_pdaf_binding_tiegcm)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error writing git version pdaf-binding-tiegcm attribute')
    istat = nf90_put_att(ncid, NF90_GLOBAL, 'git version pdaf', gitversion_pdaf)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error writing git version pdaf attribute')
    istat = nf90_put_att(ncid, NF90_GLOBAL, 'git version geodetic-fortran-utilities', gitversion_geodetic_fortran_utilities)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error writing git version geodetic-fortran-utilities attribute')
    istat = nf90_put_att(ncid, NF90_GLOBAL, 'git version esmf', gitversion_esmf)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error writing git version esmf attribute')
  end subroutine

  !> Reads a coordinate variable's values from a NetCDF file, falling back to its associated dimension if no matching variable name is found.
  subroutine read_coordinate_var(ncid, var_name, var_id, values)

    ! extern
    use netcdf

    ! tie-gcm
    use nchist_module, only: handle_ncerr

    implicit none

    ! arguments
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: var_name
    integer, intent(out) :: var_id
    real, dimension(:), allocatable, intent(inout) :: values

    ! local

    integer :: istat

    integer :: dim_id
    integer :: dim_len

    integer :: ndims
    character(len=16) :: alternative_dim_name
    integer :: alternative_dim_id(1)

    istat = nf90_inq_varid(ncid, var_name, var_id)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'read_coordinate_var')


    istat = nf90_inq_dimid( ncid, var_name, dim_id )
    if (istat /= NF90_NOERR) then
      write(*,*) 'WARNING nc file has no dimension named "', trim(var_name),&
      '" looking for alternative not matching definition of coordinate variables.'

      istat = nf90_inquire_variable(ncid, var_id, ndims=ndims)
      if(ndims==1) then
        istat = nf90_inquire_variable(ncid, var_id, dimids=alternative_dim_id)
        dim_id=alternative_dim_id(1)
        istat = nf90_inquire_dimension( ncid, dim_id, name=alternative_dim_name )
        write(*,*) 'found dimension "', trim(alternative_dim_name), '"'
      else
        write(*,*) 'no alternative was found'
        return
      end if
    end if


    istat = nf90_inquire_dimension( ncid, dim_id, len=dim_len )
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'read_coordinate_var')
    allocate(values(dim_len))
    istat = nf90_get_var(ncid, var_id, values=values)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'read_coordinate_var')

  end subroutine

  !> Reads a NetCDF time variable and converts it to seconds since the model reference epoch.
  subroutine read_time_variable_and_convert(ncid,real_time,time_var_name)

    ! extern
    use ESMF

    ! intern
    use time_module, only: convert_to_seconds_since_ref

    implicit none

    ! arguments
    integer :: ncid
    real, dimension(:), allocatable, intent(out) :: real_time
    character(len=*), optional, intent(in) :: time_var_name

    ! local
    type(ESMF_Time), dimension(:), allocatable :: time

    character(len=16) :: time_var_name_

    if(present(time_var_name)) then
      time_var_name_ = trim(time_var_name)
    else
      time_var_name_ = "time"
    end if

    call read_time_variable(ncid,time,time_var_name_)

    if(allocated(time)) then
      allocate(real_time(size(time)))
      call convert_to_seconds_since_ref(time,real_time)

      deallocate(time)
    else
      write(*,*) 'read_time_variable_and_convert failed'
    end if

  end subroutine

  !> Reads a NetCDF time variable and converts it to an array of ESMF_Time values using its units attribute.
  subroutine read_time_variable(ncid,time,time_var_name)

    ! extern
    use ESMF
    use netcdf

    ! intern
    use time_module, only: start_epoch_form_att, model_reference_epoch

    implicit none

    ! arguments
    integer :: ncid
    type(ESMF_Time), dimension(:), allocatable, intent(out) :: time
    character(len=*), optional, intent(in) :: time_var_name


    ! local
    integer :: var_id
    real, dimension(:), allocatable :: real_time

    character(len=64) :: time_unit

    type(ESMF_Time) :: nc_ref_epoch

    type(ESMF_TimeInterval) :: interval

    character(len=64) :: timestr

    real :: conversion_factor

    integer :: i, istat, rc

    character(len=16) :: time_var_name_

    if(present(time_var_name)) then
      time_var_name_ = trim(time_var_name)
    else
      time_var_name_ = "time"
    end if

    write(*,*) 'read time from nc file from variable "' // trim(time_var_name_) // '"'

    call read_coordinate_var(ncid, time_var_name_, var_id, real_time)

    if(allocated(real_time).eqv..false.)then
       write(*,*) 'read_time_variable failed'
       return
    end if

    istat =nf90_get_att(ncid, var_id, 'units', time_unit)

    call start_epoch_form_att(time_unit,nc_ref_epoch,conversion_factor)
    real_time = real_time*conversion_factor

    call ESMF_TimeGet(nc_ref_epoch, timeString=timestr)
    write(*,*) 'file reference epoch is: ', trim(timestr)

    allocate(time(size(real_time)))
    do i = 1, size(real_time)
        time(i) = nc_ref_epoch
        call ESMF_TimeIntervalSet(interval, s_r8=real_time(i), rc=rc)
        if(ESMF_LogFoundError(rc,msg="read_reg_grid_nc:ESMF_TimeIntervalSet", &
            rcToReturn=rc)) then
            call shutdown("read_reg_grid_nc:ESMF_TimeIntervalSet")
        end if

        time(i) = time(i) + interval
    end do

  end subroutine

  !> Resolves up to four nested NetCDF group names to the id of the innermost group.
  function get_group_id(ncid, grp_name_1,grp_name_2,grp_name_3,grp_name_4) result(child_id)
    implicit none

    ! arguments
    integer, intent(in) :: ncid
    character(len=*), intent(in), optional :: grp_name_1,grp_name_2,grp_name_3,grp_name_4

    ! result
    integer :: child_id

    ! local
    integer :: parent_id

    ! assign invalid id
    child_id = -1

    parent_id = ncid

    if (present(grp_name_1)) then
      call get_group_id_create_if_missing(parent_id,grp_name_1,child_id)
      if (present(grp_name_2)) then
        parent_id=child_id
        call get_group_id_create_if_missing(parent_id,grp_name_2,child_id)
        if (present(grp_name_3)) then
            parent_id=child_id
            call get_group_id_create_if_missing(parent_id,grp_name_3,child_id)
            if (present(grp_name_4)) then
              parent_id=child_id
              call get_group_id_create_if_missing(parent_id,grp_name_4,child_id)
            end if
        end if
      end if
    end if

  end function

  !> Looks up a NetCDF group by name under parent_id, creating it if it does not yet exist.
  subroutine get_group_id_create_if_missing(parent_id,name,child_id)

    implicit none

    integer, intent(in) :: parent_id
    character(len=*), intent(in) :: name
    integer, intent(out) :: child_id

    integer :: istat

    istat = nf90_inq_ncid(parent_id,name,child_id)
    if (istat /= NF90_NOERR) then
      istat = nf90_def_grp(parent_id,name,child_id)
      if (istat /= NF90_NOERR) call handle_ncerr(istat,'creating group '// name)
    end if

  end subroutine

  !> Determines the position of each name in dim_name_order among the actual dimensions of a NetCDF variable (netcdf file must be open).
  subroutine orderNetcdfVar(ncid, var_id, dim_name_order, order)

    ! tie-gcm
    use nchist_module, only: handle_ncerr

    implicit none

    ! arguments
    integer, intent(in) :: ncid
    integer, intent(in) :: var_id
    character(len=*), dimension(:), intent(in) :: dim_name_order
    integer, dimension(:), intent(inout) :: order

    ! local
    integer :: ndimsp
    logical :: found
    integer, allocatable, dimension(:) :: dimids, countp

    integer :: i,j
    integer :: istat

    character(16) :: dimname

    istat = nf90_inquire_variable(ncid, var_id, ndims=ndimsp )
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error in orderNetcdfVar: inq number of dims in var ')

    if( ndimsp .ne. size(dim_name_order) ) then
        write(*,*) 'Error in orderNetcdfVar: dimension missmatch. dim_name_order has wrong number of elements'
    end if

    allocate( dimids(ndimsp) )
    allocate( countp(ndimsp) )

    istat = nf90_inquire_variable(ncid, var_id, dimids=dimids)
    if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error in orderNetcdfVar: inq dimension ids ')

    do i = 1, ndimsp
        istat = nf90_inquire_dimension(ncid, dimids(i), dimname, countp(i))
        if (istat /= NF90_NOERR) call handle_ncerr(istat,'Error  in orderNetcdfVar: inq dimname ' // dimname)
        ! WRITE (*, '(9x, i3, 2a i9)') i, ' ' ,dimname, countp(i)

        found = .false.
        ! search dim name
        do j = 1, size(dim_name_order)
            if ( trim(dimname) == trim(dim_name_order(j)) ) then
                found = .true.
                order(i) = j
                exit
            end if
        end do

        if(found .eqv. .false.) then
            write(*,'(a,a,a)') 'Error. Could not find dimension ', dimname , 'in input dimension list'
        end if

    end do

    deallocate( dimids )
    deallocate( countp )

  end subroutine orderNetcdfVar

  !> Writes a 3D array to a new per-rank, per-ensemble-member NetCDF file for debugging.
  subroutine write_mat_to_netcdf(mat,filename)

    ! tie-gcm
    use init_module,only: istep
    use mpi_module, only: mytid

    ! intern
    use mod_parallel_pdaf, only: task_id

    implicit none

    ! arguments
    character(len=*), intent(in) :: filename
    real, intent(in) :: mat(:,:,:)

    ! local
    integer :: istat, ncid, var_id
    character(len=4) :: dim_name
    integer :: i
    integer :: dim_ids(5)
    character(len=80) :: filename_rank

    write(filename_rank,'(2a,i0.3,a1,i0.3,a1,i0.5,a3)') filename,'_',task_id,'_',mytid,'_',istep,'.nc'

    istat = nf90_create(filename_rank, NF90_NETCDF4, ncid)

    do i= 1, rank(mat)
      write(dim_name,'(a3,i1)') 'dim', i
      istat = nf90_def_dim(ncid, dim_name, size(mat,i), dim_ids(i))
    end do

    istat = nf90_def_var(ncid, 'mat', NF90_DOUBLE, &
                      dim_ids(1:rank(mat)), &
                      var_id )

    istat = nf90_put_var(ncid, var_id, mat)

    istat = nf90_close(ncid)
  end subroutine write_mat_to_netcdf

  !> Writes a 2D array to a new per-rank, per-ensemble-member NetCDF file for debugging.
  ! TODO add interface for write_mat
  subroutine write_mat_2d_to_netcdf(mat,filename)

    ! tie-gcm
    use init_module,only: istep
    use mpi_module, only: mytid

    ! intern
    use mod_parallel_pdaf, only: task_id

    implicit none

    character(len=*), intent(in) :: filename
    real, intent(in) :: mat(:,:)

    integer :: istat, ncid, var_id
    character(len=4) :: dim_name
    integer :: i
    integer :: dim_ids(5)
    character(len=80) :: filename_rank

    write(filename_rank,'(2a,i0.3,a1,i0.3,a1,i0.5,a3)') filename,'_',task_id,'_',mytid,'_',istep,'.nc'

    istat = nf90_create(filename_rank, NF90_NETCDF4, ncid)

    do i= 1, rank(mat)
      write(dim_name,'(a3,i1)') 'dim', i
      istat = nf90_def_dim(ncid, dim_name, size(mat,i), dim_ids(i))
    end do

    istat = nf90_def_var(ncid, 'mat', NF90_DOUBLE, &
                      dim_ids(1:rank(mat)), &
                      var_id )

    istat = nf90_put_var(ncid, var_id, mat)

    istat = nf90_close(ncid)
  end subroutine


end module netcdf_functionality
