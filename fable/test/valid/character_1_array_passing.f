      subroutine sub(strs1)
      character strs1(2)
      write(6, *) strs1
      end

      program prog
      character*1 strs1(2)
      strs1(1) = 'X'
      strs1(2) = 'y'
      call sub(strs1)
      end
