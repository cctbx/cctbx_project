      subroutine sub(strings)
      character*(*) strings(*)
      do i=1,2
        write(6, '(a)') strings(i)
      enddo
      end

      program prog
      character*3 strings(2)
      strings(1) = 'Abc'
      strings(2) = 'dEF'
      call sub(strings)
      end
