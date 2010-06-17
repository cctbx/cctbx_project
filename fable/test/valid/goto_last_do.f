      subroutine sub(num)
      if (num .eq. 1) then
        num = 2
        goto 10
      endif
      num = 3
   10 continue
      do i=1,num
        write(6, '(i2)') i
      enddo
      end

      program prog
      num = 1
      call sub(num)
      call sub(num)
      end
