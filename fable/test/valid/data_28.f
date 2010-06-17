      program prog
      dimension nums1(2)
      dimension nums2(2)
      data (nums1(i), i=1,2), (nums2(i), i=1,2) /1,2,3,4/
      write(6, *) (nums1(i), i=1,2), (nums2(i), i=1,2)
      end
