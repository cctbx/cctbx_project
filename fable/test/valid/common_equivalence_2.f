      program prog
      common /scr/ nc(2)
      dimension nl(3)
      equivalence (nc(2), nl(1))
      do i=1,2
        nc(i) = 20+i
      enddo
      do i=1,3
        nl(i) = 30+i
      enddo
      write(6, *) nc
      write(6, *) nl
      end
