      subroutine sub(num)
      integer base_size
      parameter(base_size=3)
      dimension nums(base_size*2)
      save nums
      do i=1,base_size*2
        if (num .eq. 1) then
          nums(i) = 73 + i
        else
          nums(i) = nums(i) + 13
        endif
      enddo
      write(6, *) nums
      end

      program prog
      do i=1,3
        call sub(i)
      enddo
      end
