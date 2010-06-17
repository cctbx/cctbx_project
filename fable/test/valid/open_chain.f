      program prog
      open(10, file='fable_tmp_661de075')
      close(10)
      open(10, file='fable_tmp_661de075',
     &  form='formatted',
     &  status='unknown')
      close(10)
      open(10, file='fable_tmp_661de075',
     &  access='sequential',
     &  form='formatted',
     &  status='new',
     &  err=10)
      goto 20
   10 write(6, '(a)') 'open err statement'
   20 continue
      close(10, status='keep', err=30)
      goto 40
   30 write(6, '(a)') 'close err statement'
   40 write(6, '(a)') 'Done.'
      end
