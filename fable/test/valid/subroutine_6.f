      subroutine sub1(cstr, str)
      character cstr*(*), str*(*)
      write(6, '(3a)') '!', str, '@'
      str = cstr
      write(6, '(3a)') '#', str, '$'
      str = cstr(2:3)
      write(6, '(3a)') '%', str, '^'
      str = cstr(3:3)
      write(6, '(3a)') '&', str, '*'
      str(2:2) = cstr(1:1)
      write(6, '(3a)') '(', str, ')'
      str(1:1) = cstr(2:4)
      write(6, '(3a)') '-', str, '+'
      str(1:2) = cstr(3:4)
      end

      subroutine sub2(str)
      character str*(*)
      str = 'xy'
      call sub1('aBcD', str)
      write(6, '(3a)') '{', str, '}'
      end

      program prog
      character str2*2
      call sub2(str2)
      end
