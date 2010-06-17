      function ifun
      ifun = 34
      end

      function jfun()
      jfun = 56
      write(6, *) jfun + 12
      end

      program prog
      write(6, *) ifun()
      write(6, *) jfun()
      end
