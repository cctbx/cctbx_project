      program prog
      equivalence(num, val)
      val = 1.2
      write(6, '(f3.1)') val
      num = 0
      if (val .eq. 1.2) stop 'equivalence failure'
      end
