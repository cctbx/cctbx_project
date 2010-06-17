      program prog
      dimension nums1(32)
      dimension nums2(33)
      data nums1 /
     &  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,
     &  17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32/
      data (nums2(i), i=1,33) /
     &  2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,
     &  18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34/
      do i=1,32
        if (nums1(i) .ne. i) stop 'FAILURE 1'
      enddo
      do i=1,33
        if (nums2(i) .ne. i+1) stop 'FAILURE 2'
      enddo
      write(6, '(a)') 'OK'
      end
