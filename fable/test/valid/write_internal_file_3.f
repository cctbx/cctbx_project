      program prog
      character buf*6
      num = -2
      buf = 'AbCdEf'
      write(buf(2:4), '(i3)') num
      write(6, '(''num = ('',a,'')'')') buf
      end
