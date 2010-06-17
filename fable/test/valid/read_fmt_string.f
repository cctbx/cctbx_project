      program prog
      character buf*2
   10 read(5, '(a)') buf
      write(6, '(3a)') '[', buf, ']'
      if (buf .ne. 'st') goto 10
      end
