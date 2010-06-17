      program prog
      do i=-1,1
        do j=-1,1
          write(6, *) i, j, sign(i, j)
        enddo
      enddo
      end
