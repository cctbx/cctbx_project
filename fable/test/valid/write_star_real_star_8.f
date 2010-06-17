      program prog
      double precision val
      do ie=-30,30
        if ((ie .ge. -4 .and. ie .le. 17) .or. mod(ie, 10) .eq. 0) then
          val = 1.2345678901234567D0*10.D0**ie
          write(6, *) val
        endif
      enddo
      end
