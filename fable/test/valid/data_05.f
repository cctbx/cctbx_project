      subroutine sub
      data num /7/
      write(6, *) num, num2
      num = num + 1
      num2 = num2 + 1
      end

      program prog
      call sub
      call sub
      call sub
      end
