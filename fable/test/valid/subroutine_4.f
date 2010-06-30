c1
      subroutine sub1(letter, num) ! c2
c3
      character letter
      if (letter(1:1) .eq. 'x') then
        num = num + 10
      endif
c4
      end ! c5

c6
      subroutine sub2(letter, num)
c7
      character letter
      call sub1(letter, num)
      if (letter(1:1) .eq. 'x') then
        num = 1 + num
      else
        num = num + 2
      endif
c8
      end

c9
      program prog
c10
      call sub2('x', num)
      write(6, '(i2)') num
      call sub2('y', num)
      write(6, '(i2)') num
c11
      end
c12
