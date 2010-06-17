      program prog
      open(10, file='fable_tmp_d41cda96', status='unknown')
      write(10, '(i1)') 5
      close(10)
      open(10, file='fable_tmp_d41cda96')
      read(10, '(i1)', end=10, err=20) num1, num2
      write(6, '(a)') 'FAILURE: end=10 expected'
      goto 20
   10 write(6, '(a)') 'OK'
   20 continue
      end
