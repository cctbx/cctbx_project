      subroutine sub1(num)
      write(6, *) 'sub1'
      if (num .eq. 1) then
        call sub2(num)
      endif
      end

      subroutine sub2(num)
      write(6, *) 'sub2'
      call sub1(num+1)
      end

      program prog
      write(6, *) 'start'
      call sub1(1)
      write(6, *) 'done'
      end
