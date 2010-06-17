      subroutine sub(num_max)
      if (num_max .gt. 41) then
        num = 41
      else
        num = num_max
      endif
      write(6, '(i2)') num
      end

      program prog
      do i=39,42
        call sub(i)
      enddo
      end
