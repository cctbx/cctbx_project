      program prog
      character ld*1(3, 2)
      character ls*3(2)
      equivalence(ld, ls)
      data ls /'Yui', 'dFg'/
      write(6, *) ld
      end
