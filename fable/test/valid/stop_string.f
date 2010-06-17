      program prog
      do i=1,2
        if (i .eq. 2) then
          stop 'Break'
        endif
        write(6, '(a,i2)') 'iteration', i
      enddo
      end
