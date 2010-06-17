      program prog
      do i=0,1
        if (i .eq. 0) then
          write(6, *) 'i is zero.'
        else
          write(6, *) 'i is not zero.'
        endif
      enddo
      end
