      subroutine sub(num, nums1, nums2)
      dimension nums1(6)
      dimension nums2(2,3)
      write(6, '(i1)') num
      do i=1,6
        write(6, '(i2)') nums1(i)
      enddo
      do i=1,2
        do j=1,3
          write(6, '(i2)') nums2(i,j)
        enddo
      enddo
      end

      program prog
      dimension nums(6)
      nums(1) = 3
      do i=2,6
        nums(i) = nums(i-1) + i
      enddo
      call sub(nums, nums, nums)
      end
