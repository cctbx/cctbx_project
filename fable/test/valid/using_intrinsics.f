      program prog
      call exercise_01
      call exercise_02
      call exercise_nint
      end

      subroutine exercise_01
      num = 84
      val = float(num)
      val = val + 0.99
      num = int(val)
      write(6, '(f4.1,1x,i2)') val, num
      end

      subroutine exercise_02
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
