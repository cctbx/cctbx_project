      program prog
      dimension vals(4)
      vals(1) = 1.2
      vals(2) = vals(1)**2
      vals(3) = 2**vals(2)
      vals(4) = vals(2)**vals(3)
      write(6, '(4f5.2)') (vals(i), i=1,4)
      end
