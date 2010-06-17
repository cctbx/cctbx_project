      program prog
      save
      common /info/ nums(2)
      equivalence (nums(1), n1), (n2, nums(2))
      n1 = 12
      n2 = 34
      write(6, *) nums
      end
