      subroutine sub1(num)
      read(5, '(i2)') num
      end

      subroutine sub2(num)
      call sub1(num)
      end

      program prog
      call sub2(num)
      write(6, '(i2)') num*2
      end
