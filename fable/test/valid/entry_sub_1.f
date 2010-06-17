      subroutine subi
      write(6, *) 'subi'
      entry subc
      write(6, *) 'subc'
      end

      program prog
      call subi
      call subc
      call subi
      call subc
      end
