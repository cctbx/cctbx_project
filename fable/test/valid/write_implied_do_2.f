      program prog
      i = 1
      write(6, '(4i2)') i, (i,i=2,3), i
      write(6, '(2i2)') (i,i=5,6)
      write(6, '(i2)') i
      write(6, '(4i2)') (i,i=0,1), (j,j=3,4)
      end
