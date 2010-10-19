      program prog
      character*1 cx
      byte bx
      equivalence(bx, cx)
      data bx /88/
      write(6, '(a)') cx
      end
