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
! TIE-GCM/ESMF time handling: model start/end/reference epoch, ISO time parsing, GPS-UTC leap-second conversion, and time search utilities.

module time_module

  ! extern
  use esmf
  use interpolation_module, only: nn_interpolator

  implicit none

  integer, parameter :: TIME_SYSTEM_UTC = 0
  integer, parameter :: TIME_SYSTEM_GPS = 1

  interface get_current_modeltime
    module procedure get_current_modeltime_as_esmftime, &
                     get_current_modeltime_as_esmfinterval, &
                     get_current_modeltime_as_secs
  end interface get_current_modeltime

  interface convert_to_seconds_since_ref
    module procedure convert_to_seconds_since_ref_array, &
                     convert_to_seconds_since_ref_scalar
  end interface convert_to_seconds_since_ref

  type(ESMF_TimeInterval) :: offset_GPS_TAI

  real, dimension(28), parameter :: leap_seconds =(/ &
    10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, &
    19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, &
    29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0 /);

  real, dimension(28), parameter :: leap_seconds_epoch =(/ &
    41317., 41499., 41683., 42048., 42413., 42778., 43144., &
    43509., 43874., 44239., 44786., 45151., 45516., 46247., 47161., &
    47892., 48257., 48804., 49169., 49534., 50083., 50630., 51179., &
    53736., 54832., 56109., 57204., 57754. /);

  real, dimension(28), parameter ::leap_seconds_epoch_tai = leap_seconds_epoch+leap_seconds/86400.;

  type(nn_interpolator) :: get_leap_second_tai

  type(ESMF_Time) :: model_start
  type(ESMF_Time) :: model_end
  type(ESMF_Time) :: model_ref_epoch

  contains

  !> Initializes leap-second data, the GPS-to-TAI offset, and the module's
  !! model_start/model_end/model_ref_epoch times, printing a summary.
  subroutine init_time_module()
    use search_module, only: NN_LEFT

    implicit none

    integer :: rc
    character(len=64), dimension(3) :: timestr

    call get_leap_second_tai%init(leap_seconds_epoch_tai,leap_seconds,NN_LEFT,assume_sorted=.true.)
    call ESMF_TimeIntervalSet(offset_GPS_TAI, s_r8=19.0, rc=rc)

    call model_start_time( model_start )
    call model_end_time( model_end )
    call model_reference_epoch(model_ref_epoch)

    call ESMF_TimeGet(model_start, timeString=timestr(1))
    call ESMF_TimeGet(model_end, timeString=timestr(2))
    call ESMF_TimeGet(model_ref_epoch, timeString=timestr(3))

    write(*,*) 'Initalize Time module.', &
               ' Model Start: ', trim(timestr(1)), &
               ' Model End: ', trim(timestr(2)), &
               ' Model Reference Epoch:', trim(timestr(3))
  end subroutine

  !> Deallocates the leap-second interpolator.
  subroutine finalize_time_module()
    call get_leap_second_tai%deallocate
  end subroutine

  !> Parses a CF-style time-units attribute (e.g. 'minutes since
  !! 2010-01-01 00:00:00') into a reference datetime and, optionally, the
  !! unit-to-seconds conversion factor.
  subroutine start_epoch_form_att( time_att_str, datetime, conversion_factor )

    implicit none

    ! args
    character(len=*), intent(in) :: time_att_str
    type(ESMF_Time), intent(out) :: datetime
    real, intent(out), optional :: conversion_factor ! factor to convert intervals to seconds

    ! local
    integer :: pos
    character(len=20) :: datetime_str

    pos = index(time_att_str, 'since')
    datetime_str = time_att_str( pos+5 : LEN_TRIM(time_att_str) )

    call iso_time_str_to_esfm_time( datetime_str, datetime )

    if(present(conversion_factor))then
      select case(trim(time_att_str(1:pos-1)))
        case('seconds')
          conversion_factor = 1.
        case('minutes')
          conversion_factor = 60.
        case('hours')
          conversion_factor = 3600.
        case('days')
          conversion_factor = 86400.
        case default
          call shutdown('start_epoch_form_att: can not convert '//time_att_str(1:pos-1))
      end select
    end if

  end subroutine start_epoch_form_att

  !> Parses date-time strings like '2010-01-01T00:00:00', '2010-01-01
  !! 00:00:00', or '2010-1-1 0:0:0.0' (fractional seconds neglected) into
  !! year, month, day, hour, minute, and second.
  subroutine parse_iso_time_str( iso_str, year, month, day, hour, minute, sec )

    implicit none

    ! args
    character(len=*), intent(inout) :: iso_str
    integer, intent(inout) :: year, month, day, hour, minute, sec

    ! local
    integer :: seperator_pos(6)
    integer :: c, i, n

    iso_str = trim( ADJUSTL(iso_str) )

    n=LEN_TRIM(iso_str)

    seperator_pos(1:5)=-1
    seperator_pos(6) = n+1

    c = 1
    do i=1,n
      select case(c)
        case(1,2)
          if(iso_str(i:i).eq.'-')then
            seperator_pos(c)=i
            c=c+1
          end if
        case(3)
          if((iso_str(i:i).eq.'T').or.(iso_str(i:i).eq.' '))then
            seperator_pos(c)=i
            c=c+1
          end if
        case(4,5)
          if(iso_str(i:i).eq.':')then
            seperator_pos(c)=i
            c=c+1
          end if
        case(6)
          if(iso_str(i:i).eq.'.')then
            seperator_pos(c)=i
            c=c+1
          end if
      end select
    end do

    if(any(seperator_pos(1:5)<1))then
      call shutdown('failed to parse date from '//iso_str)
    end if

    read (iso_str(1                 :seperator_pos(1)-1), *) year
    read (iso_str(seperator_pos(1)+1:seperator_pos(2)-1), *) month
    read (iso_str(seperator_pos(2)+1:seperator_pos(3)-1), *) day
    read (iso_str(seperator_pos(3)+1:seperator_pos(4)-1), *) hour
    read (iso_str(seperator_pos(4)+1:seperator_pos(5)-1), *) minute
    read (iso_str(seperator_pos(5)+1:seperator_pos(6)-1), *) sec

  end subroutine parse_iso_time_str

  !> Converts a date-time string like '2010-01-01T00:00:00' or
  !! '2010-01-01 00:00:00' to an ESMF_Time.
  subroutine iso_time_str_to_esfm_time( iso_str, datetime )

    implicit none

    ! args
    character(len=*), intent(inout) :: iso_str
    type(ESMF_Time), intent(out) :: datetime

    ! local
    integer :: year, month, day, hour, minute, sec
    integer :: rc
!     character(len=64) :: timestr

    call parse_iso_time_str( iso_str, year, month, day, hour, minute, sec)



    call ESMF_TimeSet(datetime,&
                      yy=int(year,kind=ESMF_KIND_I4),&
                      mm=int(month,kind=ESMF_KIND_I4),&
                      dd=int(day,kind=ESMF_KIND_I4),&
                      h=int(hour,kind=ESMF_KIND_I4),&
                      m=int(minute,kind=ESMF_KIND_I4),&
                      s=int(sec,kind=ESMF_KIND_I4),&
                      calkindflag=ESMF_CALKIND_GREGORIAN,&
                      rc=rc)
    if(ESMF_LogFoundError(rc,msg="iso_time_str_to_esfm_time:ESMF_TimeSet " // iso_str, &
        rcToReturn=rc)) then
        call shutdown('iso_time_str_to_esfm_time:ESMF_TimeSet ' // iso_str)
    end if

!     call ESMF_TimeGet(datetime, timeString=timestr)
!     write(*,'(a24,a3,i5,5i3,a3,a40)') iso_str, ' ->', year, month, day, hour, minute, sec,' ->', trim(timestr)

  end subroutine iso_time_str_to_esfm_time

  !> Returns the time of TIE-GCM's first model time step, computed from
  !! the reference epoch and pristart.
  subroutine model_start_time( start_time )

    ! tie-gcm
    use input_module, only: start_year, pristart ! ( doy, hour, minute, seconds)

    implicit none

    ! args
    type(ESMF_Time), intent(out) :: start_time

    ! local
    integer :: rc
    type(ESMF_TimeInterval) :: timeinterval

    call model_reference_epoch(start_time)

    ! substract one day. doy 1 -> no day needs to be added
    call ESMF_TimeIntervalSet(timeinterval,&
                              d=pristart(1,1)-1,&
                              h=pristart(2,1),&
                              m=pristart(3,1),&
                              s=pristart(4,1),&
                              rc=rc)
    if(ESMF_LogFoundError(rc,msg="model_start_time:ESMF_TimeIntervalSet", &
        rcToReturn=rc)) then
        call shutdown("model_start_time:ESMF_TimeIntervalSet")
    end if

    start_time = start_time + timeinterval

  end subroutine

  !> Returns the model end time, computed as the start time plus nstep
  !! TIE-GCM steps.
  subroutine model_end_time( end_time )

    ! extern
    use ESMF

    ! tie-gcm
    use hist_module, only: nstep  ! total_steps in TIE-GCM
    use input_module, only: step  ! TIE-GCM temporal resolution

    implicit none

    ! args
    type(ESMF_Time), intent(out) :: end_time

    ! local
    integer :: rc
    type(ESMF_TimeInterval) :: timeinterval

    call model_start_time(end_time)
    call ESMF_TimeIntervalSet(timeinterval, s=nstep*step, rc=rc)
    end_time = end_time + timeinterval

  end subroutine

  !> Returns the model's reference epoch, January 1 of start_year.
  subroutine model_reference_epoch(ref_time)

    ! extern
    use ESMF

    ! tie-gcm
    use input_module, only: start_year

    implicit none

    ! arguments
    type(ESMF_Time), intent(out) :: ref_time

    ! local
    integer :: rc

    call ESMF_TimeSet(ref_time,&
                      yy=start_year,&
                      mm=1,&
                      dd=1,&
                      h=0,&
                      m=0,&
                      s=0,&
                      calkindflag=ESMF_CALKIND_GREGORIAN,&
                      rc=rc)
    if(ESMF_LogFoundError(rc,msg="model_reference_epoch:ESMF_TimeSet ", &
        rcToReturn=rc)) then
        call shutdown('model_reference_epoch:ESMF_TimeSet ')
    end if
  end subroutine

  !> Returns the current model time as an absolute ESMF_Time (reference
  !! epoch plus elapsed model time).
  subroutine get_current_modeltime_as_esmftime( time )

    implicit none

    ! args
    type(ESMF_Time), intent(out) :: time

    ! local
    type(ESMF_TimeInterval) :: modeltime

    call model_reference_epoch(time)
    call get_current_modeltime(modeltime)

    time = time + modeltime

  end subroutine

  !> Returns the time elapsed since the reference epoch, from TIE-GCM's
  !! current modeltime (doy, hour, minute, second).
  subroutine get_current_modeltime_as_esmfinterval(timeinterval)

    ! tie-gcm
    use hist_module,only: modeltime ! integer(4) -> (doy,hrs,mins,secs)

    implicit none

    ! arguments
    type(ESMF_TimeInterval), intent(out) :: timeinterval

    ! local
    integer :: rc
    type(ESMF_Time) :: ref_time

    call model_reference_epoch(ref_time)

    ! substract one day. doy 1 -> no day needs to be added
    call ESMF_TimeIntervalSet(timeinterval,&
                              ref_time,&
                              d=modeltime(1)-1,&
                              h=modeltime(2),&
                              m=modeltime(3),&
                              s=modeltime(4),&
                              rc=rc)
    if(ESMF_LogFoundError(rc,msg="model_time_current_step:ESMF_TimeIntervalSet", &
        rcToReturn=rc)) then
        call shutdown("model_time_current_step:ESMF_TimeIntervalSet")
    end if
  end subroutine

  !> Returns the time elapsed since the reference epoch, in seconds.
  subroutine get_current_modeltime_as_secs(seconds)

    use ESMF

    implicit none

    ! arguments
    real(ESMF_KIND_R8), intent(out) :: seconds

    ! local
    type(ESMF_TimeInterval) :: timeinterval

    call get_current_modeltime(timeinterval)
    call ESMF_TimeIntervalGet(timeinterval,s_r8=seconds)

  end subroutine

  !> Converts an ESMF_TimeInterval to an equivalent number of TIE-GCM time
  !! steps.
  function timeinterval_to_tgcm_steps( interval ) result(steps)

    implicit none

    ! args
    type(ESMF_TimeInterval), intent(in) :: interval

    ! result
    integer :: steps

    ! local
    integer(ESMF_KIND_I4) :: seconds
    integer :: rc

    call ESMF_TimeIntervalGet(interval, s=seconds, rc=rc)
    steps = seconds_to_tgcm_steps(seconds)

  end function timeinterval_to_tgcm_steps

  !> Converts a duration in seconds to a whole number of TIE-GCM time
  !! steps, aborting if not an exact multiple.
  function seconds_to_tgcm_steps( seconds ) result(steps)

    ! tie-gcm
    use input_module, only: step  ! TIE-GCM temporal resolution in seconds

    implicit none

    ! args
    integer, intent(in) :: seconds

    ! result
    integer :: steps

    if( mod(seconds,step) .ne. 0 )then
      call shutdown('cannot convert time interval to TIE-GCM steps. Interval must be a multiple of TIE-GCM step')
    end if

    steps = seconds/step

  end function seconds_to_tgcm_steps


  !> Converts an array of ESMF_Time values to seconds elapsed since the
  !! model reference epoch.
  subroutine convert_to_seconds_since_ref_array(time,modeltime)
    ! extern
    use ESMF

    implicit none

    ! arguments
    type(ESMF_Time), dimension(:), intent(in) :: time
    real, dimension(:), intent(out) :: modeltime

    ! local
    type(ESMF_Time) :: tiegcm_ref_epoch
    type(ESMF_TimeInterval) :: interval

    character(len=64) :: timestr

    integer :: i

    call model_reference_epoch(tiegcm_ref_epoch)

    call ESMF_TimeGet(tiegcm_ref_epoch, timeString=timestr)
    write(*,*) 'model reference epoch is: ', trim(timestr)

    do i = 1, size(time)
      interval = time(i)-tiegcm_ref_epoch
      call ESMF_TimeIntervalGet(interval,s_r8=modeltime(i))
    end do
  end subroutine

  !> Converts a single ESMF_Time value to seconds elapsed since the model
  !! reference epoch.
  subroutine convert_to_seconds_since_ref_scalar(time,modeltime)
    ! extern
    use ESMF

    implicit none

    ! arguments
    type(ESMF_Time), intent(in) :: time
    real, intent(out) :: modeltime

    ! local
    type(ESMF_Time) :: tiegcm_ref_epoch
    type(ESMF_TimeInterval) :: interval

    character(len=64) :: timestr

    call model_reference_epoch(tiegcm_ref_epoch)

    call ESMF_TimeGet(tiegcm_ref_epoch, timeString=timestr)
    write(*,*) 'model reference epoch is: ', trim(timestr)

    interval = time-tiegcm_ref_epoch
    call ESMF_TimeIntervalGet(interval,s_r8=modeltime)

  end subroutine

  !> Converts an array of GPS times to UTC in place, via TAI and
  !! interpolated leap seconds; entries can be skipped with an optional
  !! mask.
  subroutine convert_GPS_to_UTC(time,mask)

    ! extern
    use ESMF

    implicit none

    ! arguments
    type(ESMF_Time), dimension(:), intent(inout) :: time
    logical, dimension(size(time)), optional, intent(in) :: mask

    ! local
    type(ESMF_TimeInterval) :: leap_seconds
    type(ESMF_Calendar) :: mjdCalendar
    type(ESMF_Time) :: time_mjd
    real :: mjd
    integer :: rc
    integer :: i
!     character(len=64) :: timestr

    mjdCalendar = ESMF_CalendarCreate(ESMF_CALKIND_MODJULIANDAY, &
                                      name="mjd", rc=rc)

    do i=1,size(time)

      if(present(mask)) then
        if (mask(i).eqv..false.) cycle
      end if

      call ESMF_TimeValidate(time(i), rc=rc)
       if(ESMF_LogFoundError(rc,msg="convert_GPS_to_UTC", &
        rcToReturn=rc)) then

        write(*,*) 'convert_GPS_to_UTC invalid time at ', i
      end if

      ! to TAI
      time(i) = time(i) + offset_GPS_TAI

      ! get TAI as mjd
      time_mjd = time(i)
      call ESMF_TimeSet(time_mjd, calendar=mjdCalendar, rc=rc)
      call ESMF_TimeGet(time_mjd, d_r8=mjd, rc=rc)

      ! TAI to UTC
      call ESMF_TimeIntervalSet(leap_seconds,&
                                s_r8=get_leap_second_tai%interpolate(real(mjd)),&
                                rc=rc)

!       write(*,*) '-----'
!       call ESMF_TimeGet(time(i), timeString=timestr)
!       write(*,*) 'TAI is: ', trim(timestr), ' mjd:', mjd
!       call ESMF_TimeIntervalGet(leap_seconds,timeString=timestr)
!       write(*,*) 'leap seconds is: ', trim(timestr), ' ', get_leap_second_tai%interpolate(mjd)

      time(i) = time(i) - leap_seconds

!       call ESMF_TimeGet(time(i), timeString=timestr)
!       write(*,*) 'UTC is: ', trim(timestr)

    end do

     call ESMF_CalendarDestroy(mjdCalendar, rc=rc)
  end subroutine

  !> Converts a modified Julian date to an ESMF_Time (Gregorian calendar).
  function mjd_to_esmf_time(mjd) result(time)

    implicit none

    ! arguments
    real, intent(in) :: mjd

    ! result
    type(ESMF_Time) :: time

    ! local
    integer :: rc
!     character(len=64) :: timestr

    call ESMF_TimeSet(                          &
        time        = time,                    &
        d_r8        = mjd,                     &
        calkindflag = ESMF_CALKIND_MODJULIANDAY, &
        rc          = rc)

    if (rc /= ESMF_SUCCESS) then
      error stop 'mjd_to_esmf_time: ESMF_TimeSet failed'
    end if

    call ESMF_TimeSet(time, calkindflag=ESMF_CALKIND_GREGORIAN, rc=rc)

!     call ESMF_TimeGet(time, timeString=timestr)
!     write(*,*) trim(timestr)

  end function

  !> Binary-searches a sorted ESMF_Time array for val, returning its index
  !! or -1 if not present.
  function binary_search_esmf_time(array,val) result(pos)

    ! extern
    use esmf

    implicit none

    ! arguments
    type(ESMF_Time), dimension(:), intent(in) :: array
    type(ESMF_Time), intent(in) :: val

    ! result
    integer :: pos

    pos = lower_bound_esmf_time(array,val)

    if(pos>size(array))then
      ! val is located after the last element
      pos= -1
    else if(array(pos)/=val)then
      ! lower_bound returns the first element that is not less than val.
      ! Hence a mismatch here means that val is not contained at all, either
      ! because it lies before the first element or in between two of them.
      pos= -1
    end if

  end function binary_search_esmf_time

  !> Returns the index of the first array element strictly greater than
  !! val (std::upper_bound analog).
  function upper_bound_esmf_time(array,val) result(first)

    ! extern
    use esmf

    implicit none

    type(ESMF_Time), dimension(:), intent(in) :: array
    type(ESMF_Time), intent(in) :: val

    integer :: first

    integer :: cnt, half, mid

    first = lbound(array,dim=1)
    cnt = size(array,dim=1)

    do while(cnt>0)
      mid=first
      half = cnt/2
      mid=mid+half
      if( (val < array(mid)) .eqv. .false. ) then
        first = mid+1
        cnt = cnt-half-1
      else
        cnt = half
      end if
    end do

  end function upper_bound_esmf_time

  !> Returns the index of the first array element not less than val
  !! (std::lower_bound analog).
  function lower_bound_esmf_time(array,val) result(first)

    ! extern
    use esmf

    implicit none

    ! arguments
    type(ESMF_Time), dimension(:), intent(in) :: array
    type(ESMF_Time), intent(in) :: val

    ! result
    integer :: first

    ! local
    integer :: cnt, half, mid

    first = lbound(array,dim=1)
    cnt = size(array,dim=1)

    do while(cnt>0)
      mid=first
      half = cnt/2
      mid=mid+half
      if( array(mid) < val ) then
        first = mid+1
        cnt = cnt-half-1
      else
        cnt = half
      end if
    end do

  end function lower_bound_esmf_time

end module time_module
