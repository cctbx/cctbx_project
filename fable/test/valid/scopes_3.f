      program prog
      parameter(n_data=10)
      REAL data(n_data)
      if (n_data .le. 100) then
        write(6, *) 'branch 1.'
      else
        d_max = 0
        DO i=1,n_data
          d = data(i)
          if (d_max .lt. d) d_max = d
        enddo
      endif
      end
