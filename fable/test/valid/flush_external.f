      program prog
      external flush
      call flush(2*5-4)
      end
      subroutine flush(i)
      write(6, *) i
      end
