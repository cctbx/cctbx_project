      subroutine show1(i)
      write(6, *) 10+i
      end

      subroutine show2(i)
      write(6, *) 20+i
      end

      subroutine show(which, i)
      external which
      call which(i)
      end

      program prog
      external show1, show2
      call show(show1, 3)
      call show(show2, 4)
      end
