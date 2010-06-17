      program prog
      character s4*4
      character s12*2
      character s34*2
      parameter(s12='xy')
      parameter(s34='ab')
      data s4(1:2) /s12/
      data s4(3:4) /s34/
      write(6, '(a)') s4
      end
