      program prog
      character s2*2
      data num, s2 /12, 'Xy'/
      write(6, '(i2, x, a)') num, s2
      end
