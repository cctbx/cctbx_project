      program prog
      do i=3,4
        if (i .lt. 4) write(6, *) 'i is less than four.'
        if (i .ge. 4) write(6, *) 'i is greater than or equal to four.'
      enddo
      end
