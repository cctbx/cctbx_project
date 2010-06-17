      program prog
      dimension nums1(2)
      dimension nums2(2)
      nums1(1) = 11
      nums1(2) = 17
      nums2(1) = 19
      nums2(2) = 23
      write(6, '(2i3)') (nums1(i) + nums2(i), i=1,2)
      end
