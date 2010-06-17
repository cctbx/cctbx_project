      program prog
      character s1*2(2), s2*3(4)
      equivalence(s1(1)(1:1), s2(2)(3:3))
      s2(1) = 'abc'
      s2(2) = 'def'
      s2(3) = 'ghi'
      s2(4) = 'jkl'
      write(6, '(a)') s1
      write(6, '(a)') s2
      s1(1) = 'PQ'
      s1(2) = 'RS'
      write(6, '(a)') s1
      write(6, '(a)') s2
      end
