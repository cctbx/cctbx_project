      program prog
      parameter(n1=110, n2=120, n3=130)
      dimension field1(n1, n2, n3)
      dimension field2(n1, n2, n3)
      dimension field3(n1, n2, n3)
      common /fields/ field1, field2, field3
      jr = 0
      call set_random(jr, field1, n1, n2, n3)
      call set_random(jr, field2, n1, n2, n3)
      call set_random(jr, field3, n1, n2, n3)
      result = 0
      do iter=1,300
        call find_max_sq(result)
      enddo
      write(6, *) result
      end

      subroutine set_random(jr, field, n1, n2, n3)
      dimension field(n1, n2, n3)
      do k=1,n3
        do j=1,n2
          do i=1,n1
            jr = mod(jr*1366+150889, 714025)
            field(i,j,k) = (mod(jr, 20000) - 10000) / 10000.0
          enddo
        enddo
      enddo
      end

      subroutine find_max_sq(result)
      parameter(n1=110, n2=120, n3=130)
      dimension field1(n1, n2, n3)
      dimension field2(n1, n2, n3)
      dimension field3(n1, n2, n3)
      common /fields/ field1, field2, field3
      do k=1,n3
        do j=1,n2
          do i=1,n1
            f = ((field1(i,j,k) - field2(i,j,k)) * field3(i,j,k))**2
            if (result .lt. f) then
              result = f
            endif
          enddo
        enddo
      enddo
      end
