      program prog
      do i=1,3
        do j=1,2
          if (i-2) 11, 12, 13
  10      if (i .eq. j) goto 14
          write(6, *) 'i != j'
          goto 14
  11      write(6, *) 'label 11'
          goto 10
  12      write(6, *) 'label 12'
          goto 10
  13      write(6, *) 'label 13'
          goto 10
  14      continue
        enddo
      enddo
      end
