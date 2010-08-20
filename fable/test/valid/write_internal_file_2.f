      program prog
      character buf*8
      dimension nums(2)
      nums(1) = -2
      nums(2) = 3
      write(buf, '(2i3)') (nums(i), i=1,2)
      write(6, '(''nums = ('',a,'')'')') buf
      call exercise_internal_file_array
      call exercise_internal_file_str_arr_ref
      end

      subroutine exercise_internal_file_array
      character*3 bufs(1)
      i = 1
      write(bufs(i), '(i1,i2)') (j,j=1,2)
      write(6, *) bufs
      end

      subroutine exercise_internal_file_str_arr_ref
      character*12 ioline(2)
      common /buffer/ ioline
      ioline(1) = 'aBcDeFgHiJkL'
      ioline(2) = 'MnOpQrStUvWx'
      write(ioline, '(i2,i1)') 12, 4
      write(6, *) ioline
      write(ioline, *) -56
      write(6, *) ioline
      end
