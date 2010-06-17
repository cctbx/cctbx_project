      subroutine sub(f)
      end

      function fun(x)
      fun = -x
      end

      program prog
      y = fun(x)
      call sub(fun)
      end
