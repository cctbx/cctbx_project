      program prog
      call exercise_file_fmt
      call exercise_internal_star
      end

      subroutine exercise_file_fmt
      open(10, file='fable_tmp_acb2b385', status='unknown')
      write(10, '(a)') 'U'
      close(10)
      open(10, file='fable_tmp_acb2b385')
      read(10, '(i1)', err=10) num
      write(6, '(a)') 'FAILURE exercise_file_fmt'
      goto 20
   10 write(6, '(a)') 'success exercise_file_fmt'
   20 continue
      end

      subroutine exercise_internal_star
      character str*5
      str = '20. *'
      read(str, *, err=10) val1, val2
      write(6, '(a)') 'FAILURE exercise_internal_star'
      goto 20
   10 write(6, '(a)') 'success exercise_internal_star'
   20 continue
      end
