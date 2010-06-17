      subroutine sub1(sub)
      implicit integer(s)
      external sub
      call sub(1)
      end

      subroutine sub2(num)
      write(6, *) num
      end

      program prog
      external sub2
      call sub1(sub2)
      end
