      program prog
      character buf*30
      character str6*6
      n_fail = 0
      open(unit=10, file='fable_tmp_33ecb352')
      write(unit=10, fmt='(i6)') -123
      close(unit=10)
      buf = '    fable_tmp_33ecb352'
      open(unit=10, file=buf, form='formatted')
      read(unit=10, fmt='(i6)') num
      close(unit=10)
      if (num .ne. -123) then
        write(6, '(a)') 'FAILURE int', num
        n_fail = n_fail + 1
      endif
      open(unit=10, file='fable_tmp_33ecb352')
      read(unit=10, fmt='(a)') str6
      close(unit=10)
      if (str6 .ne. '  -123') then
        write(6, '(a)') 'FAILURE str', str6
        n_fail = n_fail + 1
      endif
      if (n_fail .eq. 0) then
        write(6, '(a)') 'OK'
      endif
      end
