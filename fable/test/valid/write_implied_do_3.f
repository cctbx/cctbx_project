      program prog
      write(6, '(4i3)') ((i*10+j,i=1,2),j=3,4)
      write(6, '(4i4)') (((i*100+j*10+k,i=1,2),j=3,4),k=5,6)
      end
