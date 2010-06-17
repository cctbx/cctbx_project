      subroutine sub(iunit)
      write(iunit, '(a)') 'ABc'
      end

      program prog
      call sub(6)
      end
