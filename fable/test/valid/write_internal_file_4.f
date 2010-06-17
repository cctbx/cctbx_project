      program prog
      character buf*24
      buf = 'abcdefghijklmnopqrstuvwx'
      write(6, '(3a)') '[', buf, ']'
      write(buf, *) 123
      write(6, '(3a)') '[', buf, ']'
      write(buf, *) 1000000, 2000000
      write(6, '(3a)') '[', buf, ']'
      end
