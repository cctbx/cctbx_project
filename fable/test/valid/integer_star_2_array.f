      program prog
      integer*2 nums
      dimension nums(2)
      nums(1) = 9
      nums(2) = -6
      num_sum = nums(1) + nums(2)
      write(6, '(i1)') num_sum
      end
