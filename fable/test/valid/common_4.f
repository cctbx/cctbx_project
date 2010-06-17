      program prog
      common /com/ n(2)
      write(6, *) n(1), n(2)
      call sub(5)
      write(6, *) n(2), n(1)
      end

      subroutine sub(num)
      common /com/ n(2)
      n(1) = num + 1
      n(2) = num + 3
      end
