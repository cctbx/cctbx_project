      subroutine sub(num1, *, num2, *)
      if (num1+num2 .eq. 2) return 1
      if (num1 .eq. 3) return 2
      if (num1 .eq. 4) return 0
      if (num1 .eq. 5) return 3
      return
      end

      program prog
      do i=0,5
        call sub(i, *10, 1, *20)
        write(6, '(a)') 'regular'
        goto 30
   10   write(6, '(a)') 'goto 10'
        goto 30
   20   write(6, '(a)') 'goto 20'
   30   continue
      enddo
      end
