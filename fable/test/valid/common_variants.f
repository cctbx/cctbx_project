      subroutine sub1a
      common /scr/ i, j(2)
      i = 12
      j(1) = 34
      j(2) = 65
      end

      subroutine sub1b
      common /scr/ i, j, k
      write(6, '(i2,x,i2,x,i2)') i, j, k
      end

      subroutine sub2a
      common /scr/ i, x
      x = 56.78
      end

      subroutine sub2b
      common /scr/ i, x
      write(6, '(i2,x,f5.2)') i, x
      i = 91
      end

      subroutine sub3
      common /scr/ i, j
      j = 23
      end

      program prog
      call sub1a
      call sub1b
      call sub2a
      call sub2b
      call sub3
      call sub1b
      end
