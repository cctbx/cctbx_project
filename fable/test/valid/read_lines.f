      program prog
      character inp*4
      do while (.true.)
        read(5, '(a)', end=10) inp
        write(6, '(a)') inp
      enddo
   10 continue
      end
