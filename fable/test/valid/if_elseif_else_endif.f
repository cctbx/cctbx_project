      program prog
      do i=0,2
        if (i .eq. 0) then
          write(6, *) 'i is zero.'
        else if (i .eq. 1) then
          write(6, *) 'i is one.'
        else
          write(6, *) 'i is not zero or one.'
        endif
      enddo
      end
