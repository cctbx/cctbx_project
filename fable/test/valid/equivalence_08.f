      program prog
      dimension nums(2)
      equivalence (nums(1), n1), (nums(2), n2)
      n1 = 12
      n2 = 34
      write(6, *) nums
      nums(1) = 56
      nums(2) = 78
      write(6, *) n1
      write(6, *) n2
      end
