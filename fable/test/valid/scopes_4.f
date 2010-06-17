      program prog
      parameter(n_data=10)
      REAL data(n_data)
      REAL d(1)
      if (n_data .le. 100) then
        write(6, *) 'branch 1.'
      else
        d_max = 0
        DO i=1,n_data
          d(1) = data(i)
          if (d_max .lt. d(1)) d_max = d(1)
        enddo
      endif
      end
