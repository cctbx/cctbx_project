      block data globals
      common /com/ i
      data i /4/
      end
      external globals
      common /com/ i
      write(6, *) i
      end
