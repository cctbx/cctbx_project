      subroutine sub(n2, nums)
      save
      parameter(n1=2)
      integer nums(n1, n2)
      write(6, *) nums
      return
      end

      program prog
      dimension nums(2, 3)
      call sub(3, nums)
      end
