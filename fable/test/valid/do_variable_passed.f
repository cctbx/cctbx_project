      subroutine sub(iarg)
      do iarg=1,2
        write(6, *) iarg+13
      enddo
      end

      program prog
      call sub(num)
      write(6, *) num+17
      end
