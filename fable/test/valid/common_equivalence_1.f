      program prog
      common /scr/ nums(2), numx
      dimension numse(4)
      equivalence(nums, numse)
      do i=1,4
        numse(i) = 10 + i
      enddo
      write(6, *) nums(2), nums(1), nums(4), numx
      write(6, *) numse
      end
