      program prog
      parameter(n1=2)
      parameter(n2=n1*4)
      parameter(n3=(n2+n1)*5)
      parameter(n4=(n3-n2)**3)
      parameter(n5=n2*(n4+1)/2)
      parameter(n6=n1**2)
      dimension nums1(n2-5, n3-48)
      dimension nums2(n1:2, -1:1)
      dimension nums3(n6)
      write(6, *) n1
      write(6, *) n2
      write(6, *) n3
      write(6, *) n4
      write(6, *) n5
      call sub(3, 2, nums1)
      write(6, *) nums1
      call sub(1, 3, nums2)
      write(6, *) nums2
      call sub(n6, 1, nums3)
      write(6, *) nums3
      end

      subroutine sub(n1, n2, nums)
      dimension nums(n1, *)
      do i=1,n1
        do j=1,n2
          nums(i,j) = i*10+j
        enddo
      enddo
      end
