      program prog
      integer base_size
      parameter(base_size=3)
      parameter(nums_size=base_size*2)
      common /com/ nums(nums_size)
      do i=1,nums_size
        nums(i) = 13 + i * 5
      enddo
      write(6, *) nums
      end
