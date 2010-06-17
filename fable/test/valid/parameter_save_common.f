      program prog
      save
      parameter(ld=2)
      integer nums_sve(ld)
      common /cmn/ nums_cmn(ld)
      write(6, *) nums_sve
      write(6, *) nums_cmn
      end
