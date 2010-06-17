      subroutine sub1(i)
      i = 3
      write(6, *) 'sub1', i
      i = 7
      return
      end

      subroutine sub2(i, j)
      j = 4
      write(6, *) 'sub2', i, j
      i = 8
      j = 5
      return
      end

      program prog
      write(6, *) 'first line in prog.'
      call sub1(i)
      write(6, *) 'prog', i
      call sub2(i, j)
      write(6, *) 'prog', i, j
      write(6, *) 'last line in prog.'
      end
