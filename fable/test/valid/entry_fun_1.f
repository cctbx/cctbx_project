      integer function ifuni
      save num
      num = 1
      ifuni = num
      return
      entry ifunc
      num = num + num
      ifunc = num
      end

      program prog
      write(6, *) ifuni()
      write(6, *) ifunc()
      write(6, *) ifunc()
      write(6, *) ifuni()
      write(6, *) ifunc()
      write(6, *) ifunc()
      end
