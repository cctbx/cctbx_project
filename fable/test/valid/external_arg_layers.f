      program prog
      external exch_imp
      call sub(3, exch_imp)
      end

      subroutine sub(num, exch)
      implicit none
      integer num
      external exch
      call sub2(num, exch)
      end

      subroutine sub2(num, exch)
      call exch(num)
      end

      subroutine exch_imp(num)
      write(6, *) 'exch', num
      end
