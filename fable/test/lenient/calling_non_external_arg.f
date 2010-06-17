      subroutine sub1(num)
      write(6, *) num
      end

      subroutine sub2(subr)
      call subr(2)
      end

      program prog
      call sub1(1)
      call sub2(sub1)
      end
