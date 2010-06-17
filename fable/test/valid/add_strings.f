      program prog
      character a*3, b*4, c*7, d*4
      a = 'x"z'
      b = 'i\''l'
      c = a // b
      write(6, '(a)') c
      d = b // a
      write(6, '(a)') d
      d = a // b
      write(6, '(a)') d
      a = b // c
      write(6, '(a)') a
      end
