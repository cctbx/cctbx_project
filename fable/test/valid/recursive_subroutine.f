      program test
      call sub(3)
      end

      recursive subroutine sub(i)
      write(6,*) i
      if (i .ne. 0) then
        call sub(i-1)
      endif
      end
