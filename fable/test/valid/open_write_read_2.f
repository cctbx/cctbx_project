      program prog
      open(10, file='fable_tmp_7895777d',
     &  form='unformatted',
     &  status='unknown')
      write(10) -123
      close(10)
      open(10, file='fable_tmp_7895777d',
     &  form='unformatted',
     &  status='old')
      read(10) num
      close(10)
      if (num .ne. -123) then
        write(6, '(a)') 'FAILURE int', num
      else
        write(6, '(a)') 'OK'
      endif
      end
