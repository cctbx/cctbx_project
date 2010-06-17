      subroutine subi(i)
      save last_i
      last_i = i
      write(6, *) 'subi i:', i
      entry subc(i)
      last_i = last_i + 10
      i = last_i
      end

      program prog
      i = 1
      call subi(i)
      write(6, *) 'prog i:', i
      call subc(j)
      write(6, *) 'prog i:', i
      write(6, *) 'prog j:', j
      i = 2
      call subi(i)
      write(6, *) 'prog i:', i
      call subc(j)
      write(6, *) 'prog i:', i
      write(6, *) 'prog j:', j
      end
