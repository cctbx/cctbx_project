      program prog
      common /tab/ na, nb(2), nc(0:1), nd(-1:2, 3)
      dimension nums(17)
      equivalence (na, nums)
      do i=1,17
        nums(i) = 83 + i
      enddo
      write(6, *) na, nb, nc, nd
      end
