      program prog
      parameter(itwo=2)
      dimension nums1(2), nums2(itwo)
      parameter(ione=1)
      equivalence(nums1(1), nums2(ione))
      equivalence(nums1(ione), nums2(1))
      do i=1,itwo
        nums2(i) = 20 + i
      enddo
      write(6, *) nums1
      do i=1,itwo
        nums1(i) = 30 + i
      enddo
      write(6, *) nums2
      end
