      subroutine sub
      common /com/ num
      data num /7/
      write(6, *) num
      num = num + 1
      end

      program prog
      call sub
      call sub
      call sub
      end
