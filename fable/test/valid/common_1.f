      subroutine sub
      common /com/ num
      write(6, *) num
      return
      end

      program prog
      common /com/ num
      call sub
      num = 7
      call sub
      end
