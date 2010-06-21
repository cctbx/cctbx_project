      program prog
      parameter(n1=110, n2=120, n3=130)
      dimension field1(n1, n2, n3)
      dimension field2(n1, n2, n3)
      dimension field3(n1, n2, n3)
      common /fields/ field1, field2, field3
      dimension max_indices(3)
      jr = 0
      call set_random(jr, field1, n1*n2*n3)
      call set_random(jr, field2, n1*n2*n3)
      call set_random(jr, field3, n1*n2*n3)
      i = 0
      j = 0
      k = 0
      do iter=1,300
        i = mod(i + iter * 13, n1) + 1
        j = mod(j + iter * 17, n2) + 1
        k = mod(k + iter * 19, n3) + 1
        field1(i,j,k) = field2(i,j,k) * -1.3
        call find_max_sq(max_indices)
      enddo
      write(6, *) max_indices
      end

      subroutine set_random(jr, a, n)
      dimension a(*)
      do i=1,n
        jr = mod(jr*1366+150889, 714025)
        a(i) = (mod(jr, 20000) - 10000) / 10000.0
      enddo
      end

      subroutine find_max_sq(max_indices)
      dimension max_indices(3)
      parameter(n1=110, n2=120, n3=130)
      dimension field1(n1, n2, n3)
      dimension field2(n1, n2, n3)
      dimension field3(n1, n2, n3)
      common /fields/ field1, field2, field3
      max_f = 0
      do k=1,n3
        do j=1,n2
          do i=1,n1
            f = ((field1(i,j,k) - field2(i,j,k)) * field3(i,j,k))**2
            if (max_f .lt. f) then
              max_f = f
              max_indices(1) = i
              max_indices(2) = j
              max_indices(3) = k
            endif
          enddo
        enddo
      enddo
      end
