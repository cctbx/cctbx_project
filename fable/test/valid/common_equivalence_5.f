      subroutine sub1
      integer data(4)
      dimension inside(2)
      common /scr/ inside
      equivalence(data, inside)
      do i=1,4
        data(i) = 20+i
      enddo
      write(6, *) inside
      write(6, *) data
      end

      subroutine sub2
      dimension inside(3)
      common /scr/ inside
      write(6, *) inside
      write(6, *) inside(4)
      end

      program prog
      call sub1
      call sub2
      end
