      program prog
      real data(2, 3)
      do i=1,2
        do j=1,3
          data(i,j) = i+10*j
        enddo
      enddo
      do j=1,3
        do i=1,2
          write(6, *) data(i,j)
        enddo
      enddo
      end
