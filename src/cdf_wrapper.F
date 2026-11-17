#ifdef USE_CDF
! The CDF Fortran interface is not provided as a module but only via cdf.inc,
! which is imported line-by-line using include. As `cdf.inc` is formatted in
! fixed-format. This module wraps the required CDF routines so that
! free-form code can use them via 'use cdf90'.
!
! This file must be compiled with -fsecond-underscore, as libcdf exports
! its Fortran routines with two trailing underscores.
      module cdf90

      implicit none

      include 'cdf.inc'

      external CDF_get_var_num

      contains

      subroutine cdf90_open(file, id, status)
        character(len=*), intent(in) :: file
        integer(4), intent(out) :: id
        integer(4), intent(out) :: status
        call CDF_open_cdf(file, id, status)
      end subroutine

      subroutine cdf90_close(id, status)
        integer(4), intent(in) :: id
        integer(4), intent(out) :: status
        call CDF_close_cdf(id, status)
      end subroutine

      subroutine cdf90_get_zvars_maxwrittenrecnum(id, rec_num, status)
        integer(4), intent(in) :: id
        integer(4), intent(out) :: rec_num
        integer(4), intent(out) :: status
        call CDF_get_zvars_maxwrittenrecnum(id, rec_num, status)
      end subroutine

      subroutine cdf90_get_zvar_allrecords_varid(id, var_num, values,
     |                                           status)
        integer(4), intent(in) :: id
        integer(4), intent(in) :: var_num
        real(8), dimension(*), intent(out) :: values
        integer(4), intent(out) :: status
        call CDF_get_zvar_allrecords_varid(id, var_num, values, status)
      end subroutine

      subroutine cdf90_get_zvar_allrecords_varid_int1(id, var_num,
     |                                                values, status)
        integer(4), intent(in) :: id
        integer(4), intent(in) :: var_num
        integer(1), dimension(*), intent(out) :: values
        integer(4), intent(out) :: status
        call CDF_get_zvar_allrecords_varid(id, var_num, values, status)
      end subroutine

      subroutine cdf90_epoch_breakdown(epoch, year, month, day, hour,
     |                                 minute, second, msec)
        real(8), intent(in) :: epoch
        integer(4), intent(out) :: year, month, day, hour
        integer(4), intent(out) :: minute, second, msec
        call EPOCH_breakdown(epoch, year, month, day, hour,
     |                       minute, second, msec)
      end subroutine

      function cdf90_get_var_num(id, name) result(var_num)
        integer(4), intent(in) :: id
        character(len=*), intent(in) :: name

        integer(4) :: var_num

        var_num = CDF_get_var_num(id, name)
      end function

      end module
#endif
