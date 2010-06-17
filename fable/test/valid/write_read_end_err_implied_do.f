      program prog
      dimension nums(2)
      open(10, file='fable_tmp_184f95d6', status='unknown')
      write(10, '(2i3)', err=10) (i+23, i=8,9)
   10 continue
      close(10)
      open(10, file='fable_tmp_184f95d6', status='old')
      read(10, '(2i3)', end=20, err=30) (nums(i), i=1,2)
   20 continue
   30 continue
      close(10)
      if (nums(1) .ne. 31 .or. nums(2) .ne. 32) then
        write(6, '(a)') 'FAILURE'
      else
        write(6, '(a)') 'OK'
      endif
      end
