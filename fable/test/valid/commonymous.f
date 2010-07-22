      program prog
      common nums(2)
      call sub
      write(6, *) nums
      end

      subroutine sub
      common // nums(2)
      nums(1) = 1856
      nums(2) = 7893
      end
