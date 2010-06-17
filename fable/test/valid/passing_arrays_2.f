      subroutine sub(nums1, nums2)
      dimension nums1(0:1)
      dimension nums2(2:4, -1:2)
      nums1(0) = 23
      nums1(1) = 45
      nums2(4,-1) = 67
      do i=2,4
        do j=0,2
          nums2(i,j) = i*10+j
        enddo
      enddo
      end

      program prog
      dimension nums(12)
      call sub(nums, nums)
      write(6, '(6i3)') (nums(i), i=1,12)
      end
