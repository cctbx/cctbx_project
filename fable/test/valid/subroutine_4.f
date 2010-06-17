      subroutine sub1(letter, num)
      character letter
      if (letter(1:1) .eq. 'x') then
        num = num + 10
      endif
      end

      subroutine sub2(letter, num)
      character letter
      call sub1(letter, num)
      if (letter(1:1) .eq. 'x') then
        num = num + 1
      else
        num = num + 2
      endif
      end

      program prog
      call sub2('x', num)
      write(6, '(i2)') num
      call sub2('y', num)
      write(6, '(i2)') num
      end
