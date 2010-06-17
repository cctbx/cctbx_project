      program prog
      complex zero, one_two, three_four
      parameter(zero=(0.0e+0,0.0e+0))
      parameter(three=3)
      parameter(four=three+1)
      parameter(one_two=(1.0e+0,2.0e+0), three_four=(three,four))
      write(6, *) zero
      write(6, *) one_two
      write(6, *) three_four
      end
