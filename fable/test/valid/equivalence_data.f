      program prog
      common /scr/ nums(2), numx
      dimension numse(3)
      equivalence(nums, numse)
      data numse / 12, 34, 56 /
      write(6, *) numse
      write(6, *) numx, nums
      end
