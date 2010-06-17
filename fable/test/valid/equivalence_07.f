      program prog
C from section 8.2.3 of F77 standard
      CHARACTER A*4, B*4, C(2)*3
      EQUIVALENCE (A,C(1)), (B,C(2))
      a = 'abcd'
      b = 'efgh'
      write(6, '(a)') a
      write(6, '(a)') b
      write(6, '(a)') c(1)
      write(6, '(a)') c(2)
      a = 'ijkl'
      write(6, '(a)') a
      write(6, '(a)') b
      write(6, '(a)') c(1)
      write(6, '(a)') c(2)
      b = 'mnop'
      write(6, '(a)') a
      write(6, '(a)') b
      write(6, '(a)') c(1)
      write(6, '(a)') c(2)
      end
