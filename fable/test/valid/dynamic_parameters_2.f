      program prog
      parameter(nums_size=2)
      common /com/ nums(nums_size)
      call sub
      write(6, *) nums
      end

      subroutine sub
      parameter(nums_size=2)
      common /com/ nums(nums_size)
      do i=1,nums_size
        nums(i) = 13 + i * 5
      enddo
      end
