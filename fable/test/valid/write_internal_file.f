      program prog
      character buf*5
      num = -2
      write(buf, '(i3)') num
      write(6, '(''num = ('',a,'')'')') buf
      end
