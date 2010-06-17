      subroutine sub
      common /com/ num
      data num /5/
      write(6, *) num
      num = num + 1
      end

      program prog
      common /com/ num
      data num /6/
      write(6, *) num
      call sub
      write(6, *) num
      end
