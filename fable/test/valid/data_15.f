      program prog
      character*3 strings(2)
      data (strings(ix),ix=1,2) / 2*'xyz' /
      write(6, *) strings(1)
      write(6, *) strings(2)
      end
