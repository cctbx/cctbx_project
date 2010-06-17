      subroutine sub1(sub)
      call sub
      end

      subroutine sub2
      write(6, *) 'sub2'
      end

      program prog
      external sub2
      call sub1(sub2)
      end
