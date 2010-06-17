      program prog
      integer size
      parameter(size=2)
      dimension nums_local(size)
      dimension nums_save(size)
      save nums_save
      write(6, *) nums_local
      write(6, *) nums_save
      end
