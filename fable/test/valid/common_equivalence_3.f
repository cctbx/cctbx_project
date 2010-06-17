      program prog
      common /scr/ nums3
      dimension nums1(5)
      dimension nums2(3)
      dimension nums3(4)
      equivalence(nums1(2), nums2(3), nums3(4))
      do i=1,5
        nums1(i) = 20+i
      enddo
      do i=1,3
        nums2(i) = 30+i
      enddo
      do i=1,4
        nums3(i) = 40+i
      enddo
      write(6, *) nums1
      write(6, *) nums2
      write(6, *) nums3
      end
