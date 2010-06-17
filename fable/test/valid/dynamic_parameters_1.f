      program prog
      integer root_size
      parameter(root_size=3)
      parameter(nums_size=2*root_size)
      dimension nums(nums_size)
      call sub(nums)
      write(6, *) nums
      end

      subroutine sub(nums)
      integer root_size
      parameter(root_size=3)
      parameter(nums_size=2*root_size)
      dimension nums(nums_size)
      do i=1,nums_size
        nums(i) = 13+i
      enddo
      end
