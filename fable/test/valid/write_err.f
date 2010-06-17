      program prog
      open(10, file='fable_tmp_4d8a04e4', status='unknown')
      num = 7
      write(10, '(i1)', err=10) num
      write(6, '(a)') 'OK'
      goto 20
   10 write(6, '(a)') 'write err'
      goto 20
   20 continue
      end
