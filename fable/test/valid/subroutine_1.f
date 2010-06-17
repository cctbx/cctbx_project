      subroutine sub1
      write(6, *) 'output from sub1.'
      return
      end

      subroutine sub2()
      write(6, *) 'output from sub2.'
      return
      end

      program prog
      write(6, *) 'first line in prog.'
      call sub1
      call sub1()
      call sub2
      call sub2()
      write(6, *) 'last line in prog.'
      end
