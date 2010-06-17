      subroutine sub1(num)
      write(6, *) num
      return
      end

      subroutine sub2(subr)
      implicit none
      external subr
      call subr(2)
      return
      end

      program prog
      implicit none
      external sub1
      call sub1(1)
      call sub2(sub1)
      end
