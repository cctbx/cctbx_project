      function jfun()
      common/c/ i
      integer i
      jfun = 130 + i
      end

      program prog
      common/c/ i
      integer i
      integer j
      i = 1 
      write(6, *) jfun()
      i = 2
      j = jfun()
      write(6, *) 'j =', j
      i = 7
      if(jfun() == 137) then
        write(6, *) 'jfun() == 137'
      endif
      end
