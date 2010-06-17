      program prog
      character buf(2,3)*4
      i = 2
      j = 3
      buf(i,j) = '____'
      write(buf(i,j)(2:3), '(i1,a1)') i, '*'
      write(6, *) buf(i,j)
      end
