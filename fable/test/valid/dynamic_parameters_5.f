      program prog
      integer base_size
      parameter(base_size=3)
      dimension nums(base_size * 2)
      call sub(nums)
      write(6, *) nums
      end

      subroutine sub(nums)
      integer base_size
      parameter(base_size=3)
      dimension nums(base_size * 2)
      do i=1,base_size*2
        nums(i) = 61 + i * 7
      enddo
      end
