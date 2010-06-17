      program prog
      dimension nums1(0:1)
      dimension nums2(0:1, -1:2)
      dimension nums3(0:1, 3, -1:2)
      do i=0,1
        nums1(i) = 10+i
        do k=-1,2
          nums2(i,k) = 100*i+k
          do j=1,3
            nums3(i,j,k) = 200*i+30*j+k
          enddo
        enddo
      enddo
      write(6, '(2i3)') (nums1(i), i=0,1)
      write(6, '(8i4)') ((nums2(i,k), i=0,1), k=-1,2)
      write(6, '(8i5)') (((nums3(i,j,k), i=0,1), j=1,3), k=-1,2)
      end
