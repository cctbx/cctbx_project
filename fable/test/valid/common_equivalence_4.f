      program prog
      common /scr/ nums1
      dimension nums1(6)
      dimension nums2(3)
      dimension nums3(5)
      equivalence(nums1(3), nums2(2))
      equivalence(nums1(6), nums3(4))
      do i=1,6
        nums1(i) = 20+i
      enddo
      do i=1,3
        nums2(i) = 30+i
      enddo
      do i=1,5
        nums3(i) = 40+i
      enddo
      write(6, *) nums1
      write(6, *) nums2
      write(6, *) nums3
      end
