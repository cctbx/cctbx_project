      subroutine sub(nums, nums_size, nums_capacity)
      dimension nums(nums_capacity)
      do i=1,nums_size
        write(6, '(i3)') nums(i)
      enddo
      end

      program prog
      dimension nums(3)
      nums(1) = 12
      nums(2) = 34
      call sub(nums, 2, 3)
      end
