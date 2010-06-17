      subroutine sub(f)
      end

      function fun(x)
      fun = -x
      end

      program prog
      external fun
      call sub(fun)
      y = fun(x)
      end
