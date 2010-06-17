      subroutine sub
      common /com/ num, val
      write(6, *) num
      write(6, *) val
      return
      end

      program prog
      common /com/ num
      common /com/ val
      num = 3
      val = 9
      call sub
      end
