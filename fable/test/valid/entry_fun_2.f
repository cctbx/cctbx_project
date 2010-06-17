      integer function ifuni(init)
      save num
      num = init
      ifuni = num
      return
      entry ifunc(incr)
      num = num + incr
      ifunc = num
      end

      program prog
      write(6, *) ifuni(3)
      write(6, *) ifunc(2)
      write(6, *) ifunc(7)
      write(6, *) ifuni(-8)
      write(6, *) ifunc(4)
      write(6, *) ifunc(-3)
      end
