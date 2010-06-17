      subroutine sub
      integer vals
      common /com/ vals(2)
      vals(1) = vals(2) + 1
      end

      program prog
      integer vals
      common /com/ vals(2)
      vals(2) = 4
      call sub
      write(6, *) vals(1)
      end
