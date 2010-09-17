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
      call show(a // b)
      call show(b // c)
      a = 'xyz'
      c = 'abcdefg'
      c = a(1:2) // c(2:6)
      write(6, '(a)') c
      c = 'hijklmn'
      c = a(1:2) // c(2:5)
      write(6, '(a)') c
      end

      subroutine show(str)
      character str*(*)
      write(6, '(a)') str
      end
