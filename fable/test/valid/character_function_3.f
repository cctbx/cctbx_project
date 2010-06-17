      function fun(num)
      character*1 fun
      fun = 'x'
      if (num .eq. 1) fun(1:1) = 'y'
      end

      program prog
      write(6, *) fun(1)
      end
