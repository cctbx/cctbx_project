      program prog
      integer one
      parameter(one=1)
      data (num(ind), ind=1,one+one) /12,34/
      dimension num(2)
      write(6, *) num(1), num(2)
      end
