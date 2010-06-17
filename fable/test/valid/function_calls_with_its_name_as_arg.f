      function ifun(iarg)
      ifun = 12
      call sub(ifun, iarg)
      end

      subroutine sub(n1, n2)
      write(6, *) n1, n2
      end

      program prog
      write(6, *) ifun(34)
      end
