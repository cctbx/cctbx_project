      program prog
      do ix=1,3
        if (ix .eq. 1) then
          write(6, *) 'ix is one.'
        else
          ix_sum = ix_sum + ix
        endif
        write(6, *) ix_sum
      enddo
      do ix=2,3
        ix_sum_sq = ix * ix + ix_sum_sq
      enddo
      write(6, *) ix_sum_sq
      end
