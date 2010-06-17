      subroutine sub1
      common /com/ num
      data num /5/
      write(6, *) num
      num = num + 1
      end

      subroutine sub2
      common /com/ num
      data num /6/
      write(6, *) num
      num = num + 1
      end

      program prog
      call sub1
      call sub2
      end
