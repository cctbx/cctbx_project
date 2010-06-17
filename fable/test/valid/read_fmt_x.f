      program prog
      character buf2*2, buf3*3
   10 read(5, '(a2,3x,a3)') buf2, buf3
      write(6, '(5a)') '[', buf2, '][', buf3, ']'
      if (buf2 .ne. ' ') goto 10
      end
