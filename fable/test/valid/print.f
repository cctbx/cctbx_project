      program prog
      character fmt*7
      external dynfmt
      character*4 dynfmt
      print *, 12, 'Zpq'
   10 format(i2,a3)
      print 10, 34, 'Jel'
      print '(a3,i2)', 'OwM', 56
      fmt = '(a4,i1)'
      print fmt, 'TvDp', 7
      do i=1,2
        print dynfmt(i), i+8
      enddo
      print *
   20 format(3hXuW)
      print 20
      print '(i2,i3)', (i*3,i=3,4)
      end

      character*4 function dynfmt(i)
      if (i .eq. 1) then
        dynfmt = '(i1)'
      else
        dynfmt = '(i2)'
      endif
      end
