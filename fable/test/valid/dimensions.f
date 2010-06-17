      subroutine check(nums, sz)
      dimension nums(*)
      integer sz
      do i=1,sz
        if (nums(i) .ne. i) stop 'ERROR'
      enddo
      end

      subroutine sub1(nums, sz, nums_1, nums_0, nums_r)
      integer sz
      dimension nums(sz)
      dimension nums_1(*)
      dimension nums_0(0:*)
      dimension nums_r(-1:sz-2)
      ix = 1
      do i=1,sz
        nums(i) = ix
        ix = ix + 1
      enddo
      if (nums_1(13) .ne. 13) stop 'ERROR'
      if (nums_0(13) .ne. 14) stop 'ERROR'
      if (nums_r(13) .ne. 15) stop 'ERROR'
      end

      subroutine sub2(nums, sz1, sz2, nums_1, nums_0, nums_r)
      integer sz1, sz2
      dimension nums(sz1, sz2)
      dimension nums_1(sz1, *)
      dimension nums_0(0:sz1-1, 0:*)
      dimension nums_r(-1:sz1-2, 0:sz2-1)
      ix = 1
      do j=1,sz2
      do i=1,sz1
        nums(i,j) = ix
        ix = ix + 1
      enddo
      enddo
      if (sz1 .lt. 5) stop 'ERROR'
      if (sz2 .lt. 6) stop 'ERROR'
      if (nums_1(3,5) .ne. 23) stop 'ERROR'
      if (nums_0(3,5) .ne. 29) stop 'ERROR'
      if (nums_r(3,5) .ne. 30) stop 'ERROR'
      end

      subroutine sub3(nums, sz1, sz2, sz3, nums_1, nums_0, nums_r)
      integer sz1, sz2, sz3
      dimension nums(sz1, sz2, sz3)
      dimension nums_1(sz1, sz2, *)
      dimension nums_0(0:sz1-1, 0:sz2-1, 0:*)
      dimension nums_r(-1:sz1-2, 0:sz2-1, -1:sz3-2)
      ix = 1
      do k=1,sz3
      do j=1,sz2
      do i=1,sz1
        nums(i,j,k) = ix
        ix = ix + 1
      enddo
      enddo
      enddo
      if (sz1 .lt. 4) stop 'ERROR'
      if (sz2 .lt. 6) stop 'ERROR'
      if (sz3 .lt. 5) stop 'ERROR'
      if (nums_1(2,5,3) .ne. 142) stop 'ERROR'
      if (nums_0(2,5,3) .ne. 206) stop 'ERROR'
      if (nums_r(2,5,3) .ne. 263) stop 'ERROR'
      end

      subroutine sub4(nums, sz1, sz2, sz3, sz4, nums_1, nums_0, nums_r)
      integer sz1, sz2, sz3, sz4
      dimension nums(sz1, sz2, sz3, sz4)
      dimension nums_1(sz1, sz2, sz3, *)
      dimension nums_0(0:sz1-1, 0:sz2-1, 0:sz3-1, 0:*)
      dimension nums_r(-1:sz1-2, 0:sz2-1, 0:sz3-1, -1:sz4-2)
      ix = 1
      do l=1,sz4
      do k=1,sz3
      do j=1,sz2
      do i=1,sz1
        nums(i,j,k,l) = ix
        ix = ix + 1
      enddo
      enddo
      enddo
      enddo
      if (sz1 .lt. 6) stop 'ERROR'
      if (sz2 .lt. 3) stop 'ERROR'
      if (sz3 .lt. 4) stop 'ERROR'
      if (sz4 .lt. 3) stop 'ERROR'
      if (nums_1(4,2,3,1) .ne.  53) stop 'ERROR'
      if (nums_0(4,2,3,1) .ne. 187) stop 'ERROR'
      if (nums_r(4,2,3,1) .ne. 293) stop 'ERROR'
      end

      subroutine sub5(nums, sz1, sz2, sz3, sz4, sz5,
     &                      nums_1, nums_0, nums_r)
      integer sz1, sz2, sz3, sz4, sz5
      dimension nums(sz1, sz2, sz3, sz4, sz5)
      dimension nums_1(sz1, sz2, sz3, sz4, *)
      dimension nums_0(0:sz1-1, 0:sz2-1, 0:sz3-1, 0:sz4-1, 0:*)
      dimension nums_r(-1:sz1-2, 0:sz2-1, 0:sz3-1, -1:sz4-2, 0:sz5-1)
      ix = 1
      do m=1,sz5
      do l=1,sz4
      do k=1,sz3
      do j=1,sz2
      do i=1,sz1
        nums(i,j,k,l,m) = ix
        ix = ix + 1
      enddo
      enddo
      enddo
      enddo
      enddo
      if (sz1 .lt. 5) stop 'ERROR'
      if (sz2 .lt. 3) stop 'ERROR'
      if (sz3 .lt. 6) stop 'ERROR'
      if (sz4 .lt. 3) stop 'ERROR'
      if (sz5 .lt. 7) stop 'ERROR'
      if (nums_1(3,2,5,1,6) .ne. 2188) stop 'ERROR'
      if (nums_0(3,2,5,1,6) .ne. 2774) stop 'ERROR'
      if (nums_r(3,2,5,1,6) .ne. 2915) stop 'ERROR'
      end

      subroutine sub6(nums, sz1, sz2, sz3, sz4, sz5, sz6,
     &                      nums_1, nums_0, nums_r)
      integer sz1, sz2, sz3, sz4, sz5, sz6
      dimension nums(sz1, sz2, sz3, sz4, sz5, sz6)
      dimension nums_1(sz1, sz2, sz3, sz4, sz5, *)
      dimension nums_0(0:sz1-1, 0:sz2-1, 0:sz3-1, 0:sz4-1, 0:sz5-1, 0:*)
      dimension nums_r(
     &  -1:sz1-2, 0:sz2-1, -1:sz3-2, 0:sz4-1, -1:sz5-2, 0:sz6-1)
      ix = 1
      do n=1,sz6
      do m=1,sz5
      do l=1,sz4
      do k=1,sz3
      do j=1,sz2
      do i=1,sz1
        nums(i,j,k,l,m,n) = ix
        ix = ix + 1
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      if (sz1 .lt. 7) stop 'ERROR'
      if (sz2 .lt. 5) stop 'ERROR'
      if (sz3 .lt. 3) stop 'ERROR'
      if (sz4 .lt. 7) stop 'ERROR'
      if (sz5 .lt. 4) stop 'ERROR'
      if (sz6 .lt. 4) stop 'ERROR'
      if (nums_1(5,4,1,6,2,3) .ne. 11837) stop 'ERROR'
      if (nums_0(5,4,1,6,2,3) .ne. 18086) stop 'ERROR'
      if (nums_r(5,4,1,6,2,3) .ne. 19143) stop 'ERROR'
      end

      program prog
      dimension nums1(100)
      dimension nums2(5, 7)
      dimension nums3(7, 8, 5)
      dimension nums4(7, 3, 5, 4)
      dimension nums5(5, 4, 7, 3, 8)
      dimension nums6(8, 6, 3, 7, 5, 4)
      call sub1(nums1, 100, nums1, nums1, nums1)
      call check(nums1, 100)
      call sub2(nums2, 5, 7, nums2, nums2, nums2)
      call check(nums2, 5*7)
      call sub3(nums3, 7, 8, 5, nums3, nums3, nums3)
      call check(nums3, 7*8*5)
      call sub4(nums4, 7, 3, 5, 4, nums4, nums4, nums4)
      call check(nums4, 7*3*5*4)
      call sub5(nums5, 5, 4, 7, 3, 8, nums5, nums5, nums5)
      call check(nums5, 5*4*7*3*8)
      call sub6(nums6, 8, 6, 3, 7, 5, 4, nums6, nums6, nums6)
      call check(nums6, 8*6*3*7*5*4)
      write(6, '(a)') 'OK'
      end
