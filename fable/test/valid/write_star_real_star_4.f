      program prog
      do ie=-30,30
        if ((ie .ge. -4 .and. ie .le. 9) .or. mod(ie, 10) .eq. 0) then
          val = 1.2345678*10.**ie
          write(6, *) val
        endif
      enddo
      end
