      program prog
      integer size, dim
      parameter(size=3, dim=2)
      character strings*(size)(dim)
      strings(1) = 'AbC'
      strings(2) = 'dEf'
      write(6, '(i1,x,i1,x,a,x,a)') size, dim, strings
      end
