      program prog
      write(6, *) 'start'
      goto 20
   10 i = 3
      write(6, *) 'stmt 10'
      goto 30
   20 continue
      write(6, *) 'stmt 20'
      i = 2
   30 write(6, *) 'stmt 30', i
      if (i .eq. 2) goto 10
      do 40 j=1,2
        if (j .eq. 2) goto 40
        write(6, *) 'loop j is', j
   40 continue
      write(6, *) 'end'
      end
