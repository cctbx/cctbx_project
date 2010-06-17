      program prog
      common /all/ n1(2), n2(2,3)
      dimension m1(8)
      equivalence (m1, n1)
      dimension m1a(2)
      equivalence (m1a, n1(1))
      dimension m1b(2)
      equivalence (m1b(1), n1)
      dimension m1c(2)
      equivalence (m1c(1), n1(1))
      dimension m2(6)
      equivalence (m2, n2)
      dimension m2a(6)
      equivalence (m2a, n2(1,1))
      dimension m2b(6)
      equivalence (m2b(1), n2)
      dimension m2c(6)
      equivalence (m2c(1), n2(1,1))
      do i=1,8
        m1(i) = i+19
      enddo
      write(6, *) n1
      write(6, *) n2
      write(6, *) m1
      write(6, *) m1a
      write(6, *) m1b
      write(6, *) m1c
      write(6, *) m2
      write(6, *) m2a
      write(6, *) m2b
      write(6, *) m2c
      end
