      program prog
      dimension nums1(6)
      dimension nums2(2,3)
      equivalence(nums1(1), nums2(1,1))
      do i=1,6
        nums1(i) = 10+i
      enddo
      write(6, *) nums2
      do i=1,2
        do j=1,3
          nums2(i,j) = i*100+j
        enddo
      enddo
      write(6, *) nums1
      end
