      subroutine sub(iarg)
      if (iarg-2) 10,20,30
   10 write(6, *) 'smaller than two'
      return
   20 write(6, *) 'equal to two'
      return
   30 write(6, *) 'larger than two'
      end

      program prog
      do i=1,3
        call sub(i)
      enddo
      end
