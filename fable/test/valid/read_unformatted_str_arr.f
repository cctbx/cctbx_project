      program prog
      character ld*1(8)
      character ls*8
      equivalence (ls, ld)
      open(10, file='fable_tmp_8f819c56', status='unknown',
     &  form='unformatted')
      ls = 'AbcDefGh'
      write(10) ld
      close(10)
      ls = ' '
      open(11, file='fable_tmp_8f819c56', status='old',
     &  form='unformatted')
      read(11) ld
      close(11)
      if (ls .ne. 'AbcDefGh') then
        write(6, *) 'ld: [', ld, ']'
        write(6, *) 'ls: [', ls, ']'
        stop 'ERROR'
      endif
      write(6, '(a)') 'OK'
      end
