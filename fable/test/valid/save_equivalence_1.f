      subroutine sub(num)
      integer outside(3)
      dimension inside(2)
      save
      equivalence(outside(1), inside(2))
      write(6, *) inside
      write(6, *) outside
      do i=1,3
        outside(i) = num+i
      enddo
      write(6, *) inside
      write(6, *) outside
      write(6, *)
      end

      program prog
      call sub(20)
      call sub(30)
      end
