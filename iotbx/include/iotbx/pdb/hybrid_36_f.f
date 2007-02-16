C Fortran 77 port of the hy36encode() and hy36decode() functions in the
C hybrid_36.py Python prototype/reference implementation.
C See the Python script for more information.
C
C This file has no external dependencies.
C To use in your project, comment out "program exercise" at the
C bottom.
C
C This file is unrestricted Open Source (cctbx.sf.net).
C Please send corrections and enhancements to cctbx@cci.lbl.gov .
C
C See also: http://cci.lbl.gov/hybrid_36/
C
C Ralf W. Grosse-Kunstleve, Feb 2007.

C This subroutine is an implementation detail.
      subroutine encode_pure(
     &  digits,
     &  width,
     &  value,
     &  result)
      implicit none
C Input
      character digits*(*)
      integer width
      integer value
C Output
      character result*(*)
C Local
      character buf*16
      integer digits_size
      integer i, j, k, rest, val
C
      digits_size = len(digits)
      val = value
      i = 0
      j = 0
      if (val .lt. 0) then
        j = 1
        val = -val
      endif
    1 continue
        rest = val / digits_size
        k = val - rest * digits_size + 1
        i = i + 1
        buf(i:i) = digits(k:k)
        if (rest .ne. 0) then
          val = rest
          goto 1
        endif
      if (j .ne. 0) then
        i = i + 1
        buf(i:i) = '-'
      endif
      if (i .ge. width) then
        k = 1
      else
        k = width - i
        result(1:k) = ' '
        k = k + 1
      endif
    2 continue
        result(k:k) = buf(i:i)
        i = i - 1
        k = k + 1
        if (i .gt. 0) goto 2
      return
      end

C This subroutine is an implementation detail.
      subroutine decode_pure(
     &  digits_values,
     &  digits_size,
     &  s,
     &  result,
     &  diag,
     &  diag_size)
      implicit none
C Input
      integer digits_values(0:*)
      integer digits_size
      character s*(*)
C Output
      integer result
      character diag*(*)
      integer diag_size
C Local
      logical first_call
      save first_call
      data first_call /.true./
      integer iblank
      save iblank
      integer iminus
      save iminus
      integer si, dv
      logical have_minus
      logical have_non_blank
      integer value
      integer i
C
      if (first_call) then
        first_call = .false.
        iblank = ichar(' ')
        iminus = ichar('-')
      endif
      have_minus = .false.
      have_non_blank = .false.
      value = 0
      do 1 i=1,len(s)
        si = ichar(s(i:i))
        if (si .lt. 0 .or. si .gt. 127) goto 2
        if (si .eq. iblank) then
          if (.not. have_non_blank) goto 1
          value = value * digits_size
        else if (si .eq. iminus) then
          if (have_non_blank) goto 2
          have_non_blank = .true.
          have_minus = .true.
          goto 1
        else
          have_non_blank = .true.
          dv = digits_values(si)
          if (dv .lt. 0) goto 2
          value = value * digits_size
          value = value + dv
        endif
    1 continue
      if (have_minus) value = -value
      result = value
      diag_size = 0
      return
    2 diag = 'invalid number literal.'
      diag_size = 23
      return
      end

C
C hybrid-36 encoder: converts integer value to string result
C
C   width: must be 4 (e.g. for residue sequence numbers)
C               or 5 (e.g. for atom serial numbers)
C
C   value: integer value to be converted
C
C   result: string holding the result (len(result) must be >= width)
C           result(width+1:len(result)) is NOT modified in order
C           to maximize runtime performance!
C
C   diag: string holding error message, if any
C         len(diag) should be >= 80
C         diag is not modified if diag_size below is 0
C         DO NOT use if (diag .eq. ' ') to check for errors!
C         It is not supported because it would be too slow.
C
C   diag_size: length of error message, or 0 on success
C              use if (diag_size .ne. 0) to check for errors
C
      subroutine hy36encode(
     &  width,
     &  value,
     &  result,
     &  diag,
     &  diag_size)
      implicit none
C Input
      integer width
      integer value
C Output
      character result*(*)
      character diag*(*)
      integer diag_size
C Local
      character digits_upper*36
      character digits_lower*36
      data digits_upper /'0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      data digits_lower /'0123456789abcdefghijklmnopqrstuvwxyz'/
      save digits_upper
      save digits_lower
      integer i
C
      i = value
      if (width .eq. 4) then
        if (i .ge. -999) then
          if (i .lt. 10000) then
            call encode_pure(digits_upper(1:10), 4, i, result)
            diag_size = 0
            return
          endif
          i = i - 10000
          if (i .lt. 1213056) then
            i = i + 466560
            call encode_pure(digits_upper, 0, i, result)
            diag_size = 0
            return
          endif
          i = i - 1213056
          if (i .lt. 1213056) then
            i = i + 466560
            call encode_pure(digits_lower, 0, i, result)
            diag_size = 0
            return
          endif
        endif
      else if (width .eq. 5) then
        if (i .ge. -9999) then
          if (i .lt. 100000) then
            call encode_pure(digits_upper(1:10), 5, i, result)
            diag_size = 0
            return
          endif
          i = i - 100000
          if (i .lt. 43670016) then
            i = i + 16796160
            call encode_pure(digits_upper, 0, i, result)
            diag_size = 0
            return
          endif
          i = i - 43670016
          if (i .lt. 43670016) then
            i = i + 16796160
            call encode_pure(digits_lower, 0, i, result)
            diag_size = 0
            return
          endif
        endif
      else
        diag = 'unsupported width.'
        diag_size = 18
        return
      endif
      diag = 'value out of range.'
      diag_size = 19
      return
      end

C
C hybrid-36 decoder: converts string s to integer result
C
C   width: must be 4 (e.g. for residue sequence numbers)
C               or 5 (e.g. for atom serial numbers)
C
C   s: string to be converted
C        len(s) must be equal to width, or an error message is
C        returned otherwise
C
C   result: integer holding the conversion result
C
C   diag, diag_size: see hy36encode documentation above
C
      subroutine hy36decode(
     &  width,
     &  s,
     &  result,
     &  diag,
     &  diag_size)
      implicit none
C Input
      integer width
      character s*(*)
C Output
      integer result
      character diag*(*)
      integer diag_size
C Local
      character digits_upper*36
      character digits_lower*36
      data digits_upper /'0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      data digits_lower /'0123456789abcdefghijklmnopqrstuvwxyz'/
      save digits_upper
      save digits_lower
      logical first_call
      integer digits_values_upper(0:127)
      integer digits_values_lower(0:127)
      save first_call
      save digits_values_upper
      save digits_values_lower
      data first_call /.true./
      data digits_values_upper /128*-1/
      data digits_values_lower /128*-1/
      character ie_range*57
      save ie_range
      data ie_range /
     &  'internal error hy36decode: integer value out of range.'/
      integer i, di
C
      if (first_call) then
        first_call = .false.
        do 1, i=1,36
          di = ichar(digits_upper(i:i))
          if (di .lt. 0 .or. di .gt. 127) then
            diag = ie_range
            diag_size = len(ie_range)
            return
          endif
          digits_values_upper(di) = i-1
    1   continue
        do 2, i=1,36
          di = ichar(digits_lower(i:i))
          if (di .lt. 0 .or. di .gt. 127) then
            diag = ie_range
            diag_size = len(ie_range)
            return
          endif
          digits_values_lower(di) = i-1
    2   continue
      endif
      if (len(s) .eq. width) then
        di = ichar(s(1:1))
        if (di .ge. 0 .and. di .le. 127) then
          if (digits_values_upper(di) .ge. 10) then
            call decode_pure(
     &        digits_values_upper, 36, s, result, diag, diag_size)
            if (diag_size .ne. 0) return
            if (width .eq. 4) then
C                             - 10*36**(width-1) + 10**width
              result = result - 456560
              diag_size = 0
            else if (width .eq. 5) then
              result = result - 16696160
              diag_size = 0
            else
              goto 3
            endif
            return
          else if (digits_values_lower(di) .ge. 10) then
            call decode_pure(
     &        digits_values_lower, 36, s, result, diag, diag_size)
            if (diag_size .ne. 0) return
            if (width .eq. 4) then
C                             + 16*36**(width-1) + 10**width
              result = result + 756496
              diag_size = 0
            else if (width .eq. 5) then
              result = result + 26973856
              diag_size = 0
            else
              goto 3
            endif
            return
          else
            call decode_pure(
     &        digits_values_upper, 10, s, result, diag, diag_size)
            if (.not. (width .eq. 4 .or. width .eq. 5)) goto 3
            return
          endif
        endif
      endif
    3 diag = 'unsupported width.'
      diag_size = 18
      return
      end

C
C Unit tests for hy36encode and hy36decode
C
C   Can be safely ignored, but to guard yourself against
C   compiler bugs add
C
C     call tst_hybrid_36_f(.true.)
C
C   to your main program. The runtime in quick mode is in the
C   range of 0.01 seconds or less.
C
      subroutine tst_hybrid_36_f(quick)
      implicit none
C Input
      logical quick
C Local
      integer value
      character diag*80
      integer diag_size
      character s4*4
      character s5*5
      integer decoded
C
      call hy36decode(4, '    ', decoded, diag, diag_size)
      if (diag_size .ne. 0) stop 'error hy36decode empty'
      if (decoded .ne. 0) stop 'error decoded value empty'
      call hy36decode(4, '  -0', decoded, diag, diag_size)
      if (diag_size .ne. 0) stop 'error hy36decode empty'
      if (decoded .ne. 0) stop 'error decoded value -0'
      call hy36decode(4, '-999', decoded, diag, diag_size)
      if (diag_size .ne. 0) stop 'error hy36decode -999'
      if (decoded .ne. -999) stop 'error decoded value -999'
      call hy36decode(4, '9999', decoded, diag, diag_size)
      if (diag_size .ne. 0) stop 'error hy36decode 9999'
      if (decoded .ne. 9999) stop 'error decoded value 9999'
      call hy36decode(4, 'A000', decoded, diag, diag_size)
      if (diag_size .ne. 0) stop 'error hy36decode A000'
      if (decoded .ne. 10000) stop 'error decoded value A000'
      call hy36decode(4, 'ZZZZ', decoded, diag, diag_size)
      if (diag_size .ne. 0) stop 'error hy36decode ZZZZ'
      if (decoded .ne. 1223055) stop 'error decoded value ZZZZ'
      call hy36decode(4, 'a000', decoded, diag, diag_size)
      if (diag_size .ne. 0) stop 'error hy36decode a000'
      if (decoded .ne. 1223056) stop 'error decoded value a000'
      call hy36decode(4, 'zzzz', decoded, diag, diag_size)
      if (diag_size .ne. 0) stop 'error hy36decode zzzz'
      if (decoded .ne. 2436111) stop 'error decoded value zzzz'
C
      call hy36decode(5, '     ', decoded, diag, diag_size)
      if (diag_size .ne. 0) stop 'error hy36decode empty'
      if (decoded .ne. 0) stop 'error decoded value empty'
      call hy36decode(5, '   -0', decoded, diag, diag_size)
      if (diag_size .ne. 0) stop 'error hy36decode empty'
      if (decoded .ne. 0) stop 'error decoded value -0'
      call hy36decode(5, '-9999', decoded, diag, diag_size)
      if (diag_size .ne. 0) stop 'error hy36decode -9999'
      if (decoded .ne. -9999) stop 'error decoded value -9999'
      call hy36decode(5, '99999', decoded, diag, diag_size)
      if (diag_size .ne. 0) stop 'error hy36decode 99999'
      if (decoded .ne. 99999) stop 'error decoded value 99999'
      call hy36decode(5, 'A0000', decoded, diag, diag_size)
      if (diag_size .ne. 0) stop 'error hy36decode A0000'
      if (decoded .ne. 100000) stop 'error decoded value A0000'
      call hy36decode(5, 'ZZZZZ', decoded, diag, diag_size)
      if (diag_size .ne. 0) stop 'error hy36decode ZZZZZ'
      if (decoded .ne. 43770015) stop 'error decoded value ZZZZZ'
      call hy36decode(5, 'a0000', decoded, diag, diag_size)
      if (diag_size .ne. 0) stop 'error hy36decode a0000'
      if (decoded .ne. 43770016) stop 'error decoded value a0000'
      call hy36decode(5, 'zzzzz', decoded, diag, diag_size)
      if (diag_size .ne. 0) stop 'error hy36decode zzzzz'
      if (decoded .ne. 87440031) stop 'error decoded value zzzzz'
C
      if (.not. quick) then
        do 1, value=-999,2436111
          diag_size = -1
          call hy36encode(4, value, s4, diag, diag_size)
          if (diag_size .ne. 0) stop 'error hy36encode'
          decoded = -1000
          diag_size = -1
          call hy36decode(4, s4, decoded, diag, diag_size)
          if (diag_size .ne. 0) stop 'error hy36decode'
          if (decoded .ne. value) stop 'error decoded value'
    1   continue
C
        do 2, value=-9999,110000
          diag_size = -1
          call hy36encode(5, value, s5, diag, diag_size)
          if (diag_size .ne. 0) stop 'error hy36encode'
          decoded = -10000
          diag_size = -1
          call hy36decode(5, s5, decoded, diag, diag_size)
          if (diag_size .ne. 0) stop 'error hy36decode'
          if (decoded .ne. value) stop 'error decoded value'
    2   continue
      endif
C
      call hy36encode(4, -1000, s4, diag, diag_size)
      if (diag_size .eq. 0) stop 'error hy36encode range'
      if (diag .ne. 'value out of range.') stop 'error diag range'
      call hy36encode(4, 2436112, s4, diag, diag_size)
      if (diag_size .eq. 0) stop 'error hy36encode range'
      if (diag .ne. 'value out of range.') stop 'error diag range'
      call hy36decode(4, ' abc', decoded, diag, diag_size)
      if (diag_size .eq. 0) stop 'error hy36decode invalid'
      if (diag .ne. 'invalid number literal.') stop 'error diag invalid'
      call hy36decode(4, 'abc-', decoded, diag, diag_size)
      if (diag_size .eq. 0) stop 'error hy36decode invalid'
      if (diag .ne. 'invalid number literal.') stop 'error diag invalid'
      call hy36decode(4, 'A=BC', decoded, diag, diag_size)
      if (diag_size .eq. 0) stop 'error hy36decode invalid'
      if (diag .ne. 'invalid number literal.') stop 'error diag invalid'
C
      call hy36encode(5, -10000, s5, diag, diag_size)
      if (diag_size .eq. 0) stop 'error hy36encode range'
      if (diag .ne. 'value out of range.') stop 'error diag range'
      call hy36encode(5, 87440032, s5, diag, diag_size)
      if (diag_size .eq. 0) stop 'error hy36encode range'
      if (diag .ne. 'value out of range.') stop 'error diag range'
      call hy36decode(5, ' abcd', decoded, diag, diag_size)
      if (diag_size .eq. 0) stop 'error hy36decode invalid'
      if (diag .ne. 'invalid number literal.') stop 'error diag invalid'
      call hy36decode(5, 'ABCD-', decoded, diag, diag_size)
      if (diag_size .eq. 0) stop 'error hy36decode invalid'
      if (diag .ne. 'invalid number literal.') stop 'error diag invalid'
      call hy36decode(5, 'a=bcd', decoded, diag, diag_size)
      if (diag_size .eq. 0) stop 'error hy36decode invalid'
      if (diag .ne. 'invalid number literal.') stop 'error diag invalid'
C
      call hy36encode(3, 0, s5, diag, diag_size)
      if (diag_size .eq. 0) stop 'error hy36encode width'
      if (diag .ne. 'unsupported width.') stop 'error diag width'
      call hy36decode(3, '  0', decoded, diag, diag_size)
      if (diag_size .eq. 0) stop 'error hy36decode width'
      if (diag .ne. 'unsupported width.') stop 'error diag width'
      call hy36decode(3, 'A00', decoded, diag, diag_size)
      if (diag_size .eq. 0) stop 'error hy36decode width'
      if (diag .ne. 'unsupported width.') stop 'error diag width'
      call hy36decode(3, 'a00', decoded, diag, diag_size)
      if (diag_size .eq. 0) stop 'error hy36decode width'
      if (diag .ne. 'unsupported width.') stop 'error diag width'
C
      return
      end

C
C Calls unit tests above. To use this file in your project,
C comment out all lines below.
C
C To compile and run the unit tests use, e.g.:
C   f77 -o hy36 hybrid_36_f.f
C   ./hy36
C
      program exercise
      call tst_hybrid_36_f(.false.)
      end
