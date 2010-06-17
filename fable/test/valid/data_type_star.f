      program prog
      logical*1 l1
      integer*2 i2
      integer*4 i4
      integer*8 i8
      real*4 r4
      real*8 r8
      l1 = .false.
      i2 = 4
      i4 = 8
      if (i2 * 2 .eq. i4) then
        i8 = 16
        if (i4 * 2 .eq. i8) then
          write(6, '(a)') 'OK integers'
          l1 = .true.
        endif
      endif
      if (.not. l1) then
        write(6, '(a)') 'FAILURE integers'
      endif
      r4 = 3.14
      r8 = 6.28
      if (abs(r4*2 - r8) .lt. 1.e-5) then
        write(6, '(a)') 'OK reals'
      else
        write(6, '(a)') 'FAILURE reals'
      endif
      end
