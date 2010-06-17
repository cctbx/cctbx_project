      program prog
      open(1, file='fable_tmp_103fa7e8', status='unknown')
      write(1, '(i3)') 12+34
      rewind(1)
      read(1, '(i3)') num
      write(6, *) num
      rewind(unit=1, err=10)
      goto 20
   10 write(6, *) 'rewind FAILURE'
   20 continue
      write(1, '(i3)') 56+78
      rewind 1
      read(1, '(i3)') num
      write(6, *) num
      rewind 1
      endfile(1)
      rewind 1
      read(1, '(i3)', end=30) num
      write(6, *) 'endfile FAILURE'
   30 continue
      end
