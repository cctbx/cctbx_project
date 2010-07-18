      character s*5
      s = 'Abc'
      write(*, '(2a)') s, 'P'
      write(*, '(2a)') s(1:4), 'Q'
      write(*, '(2a)') s(1:3), 'R'
      write(*, '(2a)') s(1:2), 'S'
      write(*, '(2a)') s(1:len_trim(s)), 'T'
      write(*, '(2a)') s(1:len_trim(s(1:4))), 'U'
      write(*, '(2a)') s(1:len_trim(s(1:3))), 'V'
      write(*, '(2a)') s(1:len_trim(s(1:2))), 'W'
      s = ' '
      write(*, '(2a)') s(1:len_trim(s)), 'X'
      end
