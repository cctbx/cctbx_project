      program prog
      open(10, file='fable_tmp_95c8b3ac', status='unknown')
      write(10, '(i1)') 5
      close(10)
      open(10, file='fable_tmp_95c8b3ac')
      read(10, '(i1)', end=10) num1, num2
      write(6, '(a)') 'FAILURE: end=10 expected'
      goto 20
   10 write(6, '(a)') 'OK'
   20 continue
      end
