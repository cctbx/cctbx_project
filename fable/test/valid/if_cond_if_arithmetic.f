      program prog
      do i=1,2
        do j=1,3
          if (i .eq. 1) if (j-2) 10,20,30
          write(6, *) i, j
          goto 40
10        write(6, *) 'l10'
          goto 40
20        write(6, *) 'l20'
          goto 40
30        write(6, *) 'l30'
40        continue
        enddo
      enddo
      end
