      program prog
      character buf*8
      dimension nums(2)
      nums(1) = -2
      nums(2) = 3
      write(buf, '(2i3)') (nums(i), i=1,2)
      write(6, '(''nums = ('',a,'')'')') buf
      end
