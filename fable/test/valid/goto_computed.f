      program prog
      do j=1,2
        do i=1,2
          if (j .eq. 1) then
            goto (10, 20) i
          else
            goto (10, 20), i
          endif
10        write(6, *) 'statement 10', j
          goto 30
20        write(6, *) 'statement 20', j
30        continue
        enddo
      enddo
      end
