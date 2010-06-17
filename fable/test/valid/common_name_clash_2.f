      subroutine sub1init
      common /cmn1/ num1, num2
      num1 = 12
      num2 = 34
      end

      subroutine sub2init
      common /cmn2/ num2(2), num3
      num2(1) = 56
      num2(2) = 78
      num3 = 90
      end

      subroutine sub1show
      common /cmn1/ num1, num2
      write(6, *) num1, num2
      end

      subroutine sub2show
      common /cmn2/ num2(2), num3
      do i=1,2
        write(6, *) i, num2, num3
      enddo
      do i=3,4
        write(6, *) i, num2, num3
      enddo
      end

      program prog
      call sub1init
      call sub2init
      call sub1show
      call sub2show
      end
