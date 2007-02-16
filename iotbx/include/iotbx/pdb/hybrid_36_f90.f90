! Fortran 90 port of the hy36encode() and hy36decode() functions in the
! hybrid_36.py Python prototype/reference implementation.
! See the Python script for more information.
!
! This file has no external dependencies.
! To use in your project, comment out "program exercise" at the
! bottom.
!
! This file is unrestricted Open Source (cctbx.sf.net).
! Please send corrections and enhancements to cctbx@cci.lbl.gov .
!
! See also: http://cci.lbl.gov/hybrid_36/
!
! Joe Krahn, Feb 2007.
!
! hybrid-36 encoder: converts integer value to string result
!
!   width: must be 4 (e.g. for residue sequence numbers)
!               or 5 (e.g. for atom serial numbers)
!
!   value: integer value to be converted
!
!   errmsg: string holding error message, if any
!         len(errmsg) should be >= 80, or the message may be truncated.
!         errmsg is not modified unless errmsg_len is returned non-zero
!         DO NOT use "if (errmsg == ' ')" to check for errors!
!         It will work if errmsg is blanked before calling, but it is inefficient.
!
!   errmsg_len: length of error message, or 0 on success
!              use "if (errmsg_len /= 0)" to check for errors

! The hy36encode() function requires an interface because it uses an
! automatic result length defined by the caller.
! You can use an interface via a module, or put the interface in an INCLUDE file.

module hybrid36_interface_module
  interface function
    function hy36encode(value, width, errmsg, errmsg_len) result(encoded)
      implicit none
      integer, intent(in) :: value, width
      character(len=*), intent(inout) :: errmsg
      integer, intent(out) :: errmsg_len
      character(len=width) :: encoded
    end function hy36encode
    function hy36decode(string, errmsg, errmsg_len) result(decoded)
      implicit none
      character(len=*), intent(in) :: string
      character(len=*), intent(inout) :: errmsg
      integer, intent(out) :: errmsg_len
     integer :: decoded
    end function hy36decode
  end interface
end module hybrid36_interface_module

function hy36encode(value, width, errmsg, errmsg_len) result(encoded)
  implicit none
! Input
  integer, intent(in) :: width
  integer, intent(in) :: value
! Output
  character(len=width) :: encoded
  character(len=*), intent(inout) :: errmsg
  integer, intent(out) :: errmsg_len
! Local
  character(len=*), parameter :: digits_upper = &
                     '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  character(len=*), parameter :: digits_lower = &
                     '0123456789abcdefghijklmnopqrstuvwxyz'
  character(len=*), parameter :: range_error='value out of range.'
  character(len=*), parameter :: width_error='unsupported width.'

  errmsg_len = 0
  select case(width)

  case(4)
    select case(value)
    case(-999:9999)
      call encode(digits_upper(1:10), value, encoded)
    case(10000:1223055)
      call encode(digits_upper, value+456560, encoded)
    case(1223056:2436111)
      call encode(digits_lower, value-756496, encoded)
    case default
      errmsg=range_error
      errmsg_len=len(range_error)
    end select

  case(5)
    select case(value)
    case(-9999:99999)
      call encode(digits_upper(1:10), value, encoded)
    case(100000:43770015)
      call encode(digits_upper, value+16696160, encoded)
    case(43770016:87440031)
      call encode(digits_lower, value-26973856, encoded)
    case default
      errmsg=range_error
      errmsg_len=len(range_error)
    end select
  case default
    errmsg=width_error
    errmsg_len=len(width_error)
  end select
contains
! Note ranges are checked above, so we do not have to check for
! overflowing the bounds of the result string.
  subroutine encode(digits, value, result)
    implicit none
! Input
    character(len=*), intent(in) :: digits
    integer, intent(in) :: value
! Output
    character(len=*), intent(out) :: result
! Local
    integer :: i, k, n, rest, val
    logical :: negative

    n=len(digits)
    val = abs(value)
    negative = (value < 0)
    i = len(result)
    do while (val/=0)
      rest = val / n
      k = val - rest * n + 1
      result(i:i) = digits(k:k)
      i = i - 1
      val = rest
    end do
    if (negative) then
      result(i:i) = '-'
      i = i - 1
    end if
    result(1:i)=' '
    return
  end subroutine encode
end function hy36encode


! hybrid-36 decoder: converts string s to integer result
!
!   width: must be 4 (e.g. for residue sequence numbers)
!               or 5 (e.g. for atom serial numbers)
!
!   s: string to be converted
!
!   s_size: length of string to be converted
!           s_size must be <= len(s)
!           s(s_size+1:len(s)) is ignored
!
!   result: integer holding the conversion result
!
!   errmsg, errmsg_len: see hy36encode documentation above

function hy36decode(string, errmsg, errmsg_len) result(decoded)
  implicit none
! Input
  character(len=*), intent(in) :: string
! Output
  character(len=*), intent(inout) :: errmsg
  integer, intent(out) :: errmsg_len
  integer :: decoded
! Local
  character(len=*), parameter :: digits_upper = &
                     '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  character(len=*), parameter :: digits_lower = &
                     '0123456789abcdefghijklmnopqrstuvwxyz'
  logical, save :: first_call=.true.
  integer, save :: digits_values_upper(0:127) = -1
  integer, save :: digits_values_lower(0:127) = -1
  character(len=*), parameter :: ie_range = &
      'internal error hy36decode: integer value out of range.'
  character(len=*), parameter :: range_error='value out of range.'
  character(len=*), parameter :: width_error='unsupported width.'
  character(len=*), parameter :: bad_char_error= &
                                 'unsupported character in encoded string'
  integer :: i, di, width

  width=len(string)
  if (width/=4 .and. width/=5) then
    errmsg = 'unsupported width.'
    errmsg_len = 18
    return
  end if
  errmsg_len=0

  if (first_call) then
    first_call = .false.
    do i=1,36
      di = ichar(digits_upper(i:i))
      if (di < 0 .or. di > 127) then
        errmsg = ie_range
        errmsg_len = len(ie_range)
        return
      end if
      digits_values_upper(di) = i-1

      di = ichar(digits_lower(i:i))
      if (di < 0 .or. di > 127) then
        errmsg = ie_range
        errmsg_len = len(ie_range)
        return
      end if
      digits_values_lower(di) = i-1
    end do
  end if ! (first_call)

  di = ichar(string(1:1))
  if (di<0 .or. di >127) then
    errmsg = bad_char_error
    errmsg_len = len(bad_char_error)
    return
  end if

  if (digits_values_upper(di) >= 10) then
    call decode(digits_values_upper, 36, string, decoded)
    if (errmsg_len == 0) decoded = decoded - merge(456560,16696160,width==4)
  else if (digits_values_lower(di) >= 10) then
    call decode(digits_values_lower, 36, string, decoded)
    if (errmsg_len == 0) decoded = decoded + merge(756496,26973856,width==4)
  else
    call decode(digits_values_upper, 10, string, decoded)
  end if
  return
contains
  subroutine decode(digits_values, digits_base, string, result)
    implicit none
! Input
    integer, intent(in) :: digits_values(0:127)
    integer, intent(in) :: digits_base
    character(len=*), intent(in) :: string
! Output
    integer, intent(out) :: result
! Used by host association
    !character(len=*) :: errmsg
    !integer :: errmsg_len
! Local
    logical negative
    integer :: i, si, dv
    character(len=*), parameter :: invalid_literal_error = 'invalid number literal.'
! begin
    negative=.false.
    i=1
  ! First remove leading blanks, checking for bad characters,
  ! and check if the first non-blank is '-'
    do while (i<=len(string))
      select case(ichar(string(i:i)))
      case (:-1,128:);
        goto 1
      case (ichar(' '))
        continue
      case (ichar('-'))
        negative=.true.
        i=i+1
        exit
      case default
        exit
      end select
      i=i+1
    end do

    result = 0
    do while (i<=len(string))
      si=ichar(string(i:i))
      select case(si)
      case (:-1,128:)
        goto 1
      case (ichar(' '))
        result = result * digits_base
      case default
        dv = digits_values(si)
        if (dv < 0 .or. dv >= digits_base) goto 1  ! BUG FIX
        result = result * digits_base + dv
      end select
      i=i+1
    end do

    if (negative) result = -result
    errmsg_len = 0
    return

  1 continue
    errmsg = invalid_literal_error
    errmsg_len = len(invalid_literal_error)
    return
  end subroutine decode
end function hy36decode

!
! Unit tests for hy36encode and hy36decode
!
!   Can be safely ignored, but to guard yourself against
!   compiler bugs add
!
!     call tst_hybrid_36_f(.true.)
!
!   to your main program. The runtime in quick mode is in the
!   range of 0.01 seconds or less.
!
subroutine tst_hybrid_36_f(quick)
  use hybrid36_interface_module
  implicit none
! Input
  logical, intent(in) :: quick
! Local
  integer, parameter :: stderr=0
  integer :: value
  character(len=80) :: errmsg
  integer :: errmsg_len
  character(len=4) :: s4
  character(len=5) :: s5
  integer :: decoded

  call assert_both('    ',0)
  call assert_decoded('  -0',0)
  call assert_both('-999',-999)
  call assert_both('9999',9999)
  call assert_both('A000',10000)
  call assert_both('ZZZZ',1223055)
  call assert_both('a000',1223056)
  call assert_both('zzzz',2436111)

  call assert_both('     ',0)
  call assert_decoded('   -0',0)
  call assert_both('-9999',-9999)
  call assert_both('99999',99999)
  call assert_both('A0000',100000)
  call assert_both('ZZZZZ',43770015)
  call assert_both('a0000',43770016)
  call assert_both('zzzzz',87440031)

  call assert_encoded(-1000,s4,'value out of range.')
  call assert_encoded(2436112,s4,'value out of range.')
  call assert_decoded(' abc',assert_errmsg='invalid number literal.')
  call assert_decoded('abc-',assert_errmsg='invalid number literal.')
  call assert_decoded('A=BC',assert_errmsg='invalid number literal.')

  call assert_encoded(-10000,s5,'value out of range.')
  call assert_encoded(87440032,s5,'value out of range.')
  call assert_decoded(' abcd',assert_errmsg='invalid number literal.')
  call assert_decoded('ABCD-',assert_errmsg='invalid number literal.')
  call assert_decoded('a=bcd',assert_errmsg='invalid number literal.')

  call assert_encoded(0,'123','unsupported width.')
  call assert_decoded('  0',assert_errmsg='unsupported width.')
  call assert_decoded('A00',assert_errmsg='unsupported width.')
  call assert_decoded('a00',assert_errmsg='unsupported width.')

  if (quick) return

  do value=-999,2436111
    errmsg_len = -1
    s4 = hy36encode(value, 4, errmsg, errmsg_len)
    if (errmsg_len /= 0) then
      write(stderr,'(2A)') 'ERROR in hy36encode: error=',errmsg(:errmsg_len)
      stop 31
    end if
    decoded = -1000
    errmsg_len = -1
    decoded=hy36decode(s4, errmsg, errmsg_len)
    if (errmsg_len /= 0) then
      write(stderr,'(2A)') 'ERROR in hy36decode: error=',errmsg(:errmsg_len)
      stop 32
    end if
    if (value/=decoded) then
      write(stderr,*) 'ERROR in double conversion: ',value,' => "',s4,'" => ',decoded
      stop 33
    end if
  end do

  do value=-9999,110000
    errmsg_len = -1
    s5 = hy36encode(value, 5, errmsg, errmsg_len)
    if (errmsg_len /= 0) then
      write(stderr,'(2A)') 'ERROR in hy36encode: error=',errmsg(:errmsg_len)
      stop 41
    end if
    decoded = -10000
    errmsg_len = -1
    decoded = hy36decode(s5, errmsg, errmsg_len)
    if (errmsg_len /= 0) then
      write(stderr,'(2A)') 'ERROR in hy36decode: error=',errmsg(:errmsg_len)
      stop 42
    end if
    if (value/=decoded) then
      write(stderr,*) 'ERROR in double conversion: ',value,' => "',s5,'" => ',decoded
      stop 43
    end if
  end do
  return
contains
  subroutine assert_both(encoded,decoded)
    implicit none
    character(len=*), intent(in) :: encoded
    integer, intent(in) :: decoded
    call assert_encoded(decoded,encoded)
    call assert_decoded(encoded,decoded)
  end subroutine assert_both
  subroutine assert_decoded(value,assert_result,assert_errmsg)
    implicit none
    character(len=*), intent(in) :: value
    integer, intent(in), optional :: assert_result
    character(len=*), intent(in), optional :: assert_errmsg
  ! local
    integer :: result
    character(len=80) :: errmsg
    integer :: errmsg_len
  ! begin
    result=hy36decode(value, errmsg, errmsg_len)
    if (present(assert_errmsg)) then
      if (errmsg_len /= len(assert_errmsg) &
          .or. errmsg /= assert_errmsg) then
        write(stderr,*) 'ERROR in hy36decode: incorrect error message for "',value,'"'
        write(stderr,*) 'expected error "',assert_errmsg,'", got "', &
                        errmsg(:errmsg_len),'", with result=',result
        stop 11
      end if
    else ! do not check result if an error was expected
      if (errmsg_len/=0) then
        write(stderr,*) 'ERROR in hy36decode: error decoding "',value,'"'
        write(stderr,*) '      error message = "',errmsg(:errmsg_len),'"'
        stop 12
      end if
      if (present(assert_result)) then
        if (result /= assert_result) then
          write(stderr,*) 'ERROR in hy36decode: incorrect result for "',value,'"'
          write(stderr,*) 'expected result ',assert_result,', got ',result
          stop 13
        end if
      end if
    end if
  end subroutine assert_decoded
! assert_result is required, even for error assertions, because it defines the width.
  subroutine assert_encoded(value,assert_result,assert_errmsg)
    implicit none
    integer, intent(in) :: value
    character(len=*), intent(in) :: assert_result
    character(len=*), intent(in), optional :: assert_errmsg
  ! local
    character(len=len(assert_result)) :: result
    character(len=80) :: errmsg
    integer :: errmsg_len
  ! begin
    result=hy36encode(value, len(assert_result), errmsg, errmsg_len)
    if (present(assert_errmsg)) then
      if (errmsg_len /= len(assert_errmsg) &
          .or. errmsg /= assert_errmsg) then
        write(stderr,*) 'ERROR in hy36encode: incorrect error message for value=',value
        write(stderr,*) 'expected error "',assert_errmsg,'", got "', &
                        errmsg(:errmsg_len),'", with result="',result,'"'
        stop 21
      end if
    else ! do not check result if an error was expected
      if (errmsg_len/=0) then
        write(stderr,*) 'ERROR in hy36encode: error encoding value=',value
        write(stderr,*) '      error message = "',errmsg(:errmsg_len),'"'
        stop 22
      else if (result /= assert_result) then
        write(stderr,*) 'ERROR in hy36encode: incorrect result for value=',value
        write(stderr,*) 'expected result "',assert_result,'", got "',result,'"'
        stop 23
      end if
    end if
  end subroutine assert_encoded
end subroutine tst_hybrid_36_f

!
! Calls unit tests above. To use this file in your project,
! comment out all lines below.
!
! To compile and run the unit tests use, e.g.:
!   f90 -o hy36 hybrid_36_f.f90
!   ./hy36
!
program exercise
  call tst_hybrid_36_f(.false.)
end program exercise
