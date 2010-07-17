      program prog
      i = 123
      do while (i .lt. 169)
        write(6, *) i
        i = i + 45
      enddo
      do 10 while (i .lt. 281)
        write(6, *) i
        i = i + 67
   10 enddo
      end
