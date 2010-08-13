      program prog
      call exercise_c_keywords
      call exercise_char_ichar
      call exercise_nint
      call exercise_double
      call exercise_string
      call exercise_index
      end

      subroutine exercise_c_keywords
      num = 84
      val = float(num)
      val = val + 0.99
      num = int(val)
      write(6, '(f4.1,1x,i2)') val, num
      end

      subroutine exercise_char_ichar
      character c*1
      c = 'R'
      i = ichar(c)
      write(6, '(i3)') i
      c = 'X'
      write(6, '(a)') c
      c = char(i)
      write(6, '(a)') c
      end

      subroutine exercise_nint
      real rv(8)
      double precision dv(8)
      data rv /0.49, 0.51, -0.49, -0.51, 3.49, 3.51, -2.49, -2.51/
      data dv /1.49, 2.51, -3.49, -4.51, 5.49, 6.51, -7.49, -8.51/
      write(6, '(8i3)') (nint(rv(i)), i=1,8)
      write(6, '(8i3)') (nint(dv(i)), i=1,8)
      end

      subroutine exercise_double
      write(6, '(f6.2)') dacos(0.5d0)
      write(6, '(f6.2)') dasin(0.5d0)
      write(6, '(f6.2)') datan2(0.1d0, 0.2d0)
      write(6, '(f6.2)') dtan(0.5d0)
      end

      subroutine exercise_string
      character s4*4
      s4 = 'Xy'
      write(6, '(i1)') lnblnk(s4)
      end

      subroutine exercise_index
      character digits*10
      data digits /'0123456789'/
      write(6, *) index(digits, 'x')
      write(6, *) index(digits, '0')
      write(6, *) index(digits, '5')
      write(6, *) index(digits, '9')
      write(6, *) index(digits, '24')
      write(6, *) index(digits, '34')
      end
