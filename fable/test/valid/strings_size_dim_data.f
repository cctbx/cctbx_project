      program prog
      integer base, size, dim
      parameter(base=4)
      parameter(dim=base-2, size=base-1)
      character strings*(size)(dim)
      data strings /'aBc', 'DeF'/
      write(6, '(i1,x,i1,x,a,x,a)') size, dim, strings
      end
