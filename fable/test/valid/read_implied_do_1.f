      subroutine sub(nums)
      dimension nums(*)
      read(5, '(2i3)') (nums(i), i=1,2)
      end

      program prog
      dimension nums(2)
      call sub(nums)
      write(6, '(2i5)') nums(1)*3, nums(2)*4
      end
