      function fun(num)
      character fun*1
      fun = 'x'
      if (num .eq. 1) fun = 'y'
      end

      program prog
      external fun
      character fun*1
      do num=1,2
        write(6, *) fun(num)
      enddo
      end
