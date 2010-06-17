      subroutine sub
      integer num
      save num
      write(6, *) num
      num = num + 1
      return
      end

      program prog
      call sub
      call sub
      end
