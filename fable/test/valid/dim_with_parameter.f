      program prog
      parameter(num=2)
      dimension nums1(num)
      dimension nums2(num)
      nums1(1) = 12
      nums1(2) = 34
      nums2(1) = 56
      nums2(2) = 78
      write(6, '(4i3)') nums1(1), nums1(2), nums2(1), nums2(2)
      end
