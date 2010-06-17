      integer function fun(num)
      write(6, *) num
      fun = num + 7
      end

      program prog
      integer fun
      do i=1,3
        num = fun(num)
      enddo
      end
