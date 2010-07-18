      do i=1,5
        if (i .eq. 2) cycle
        write(6, '(i1)') i
      enddo
      do i=1,5
        if (i .eq. 4) exit
        write(6, '(i1)') i
      enddo
      do i=1,5
        if (i .eq. 2) cycle
        if (i .eq. 4) exit
        write(6, '(i1)') i
      enddo
      end
