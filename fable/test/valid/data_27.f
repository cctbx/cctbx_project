      program prog
      character s2s*2(2)
      data s2s(1)(1:1) /'A'/
      data s2s(2)(2:2) /'b'/
      data s2s(1)(2:2) /'C'/
      data s2s(2)(1:1) /'d'/
      write(6, '(5a)') '[', s2s(1), '][', s2s(2), ']'
      end
