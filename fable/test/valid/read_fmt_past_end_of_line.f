      program prog
      character title*2
      open(1, file='fable_tmp_77e68e5d', status='unknown')
      write(1, '(1X,i4)') 2
      write(1, '(1x,a)') 'L1'
      write(1, '(1x,a)') 'L2'
      write(1, '(1x,a)') 'L3'
      close(1)
      open(1, file='fable_tmp_77e68e5d', status='old')
      read(1, '(1x,i5)') num
      write(6, *) num
      do i=1,num
        read(1, '(1X,A2)') title
        write(6, *) '[', title, ']'
      enddo
      rewind(1)
      read(1, '(1x,f5.0)') val
      write(6, *) val
      do i=1,num
        read(1, '(1X,A2)') title
        write(6, *) '[', title, ']'
      enddo
      rewind(1)
      num = 12
      missing = 34
      read(1, '(1x,2i4)') num, missing
      write(6, *) num, missing
      close(1)
      end
