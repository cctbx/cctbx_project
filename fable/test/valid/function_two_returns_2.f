      integer function fun(i)
      if (i .lt. 3) then
        fun = i*4
        return
      endif
      fun = -i
      end

      program prog
      integer fun
      write(6, *) fun(2)
      write(6, *) fun(3)
      end
