      subroutine sub(i, j)
      if (i .lt. 2) then
        if (j .lt. 2) then
          write(6, '(a)') 'A'
          goto 10
        endif
        write(6, '(a)') 'B'
  10    continue
      else
        write(6, '(a)') 'C'
      endif
      end

      program prog
      do i=1,2
        do j = 1,2
          call sub(i, j)
        enddo
      enddo
      end
