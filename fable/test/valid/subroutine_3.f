      subroutine sub1(num)
      num = 9
      return
      end

      subroutine sub2(nums, nums_size)
      integer nums(*)
      do i=1,nums_size
        nums(i) = i * 10
      enddo
      return
      end

      subroutine sub3(nums, nums_size)
      dimension nums(*)
      do i=1,nums_size
        nums(i) = i * 20
      enddo
      return
      end

      program prog
      integer nums(2)
      call sub1(num)
      write(6, *) 'num after sub1:', num
      n = 1
      call sub2(num, n)
      write(6, *) 'num after sub2', num
      call sub1(nums)
      write(6, *) 'nums after sub1:', nums(1), nums(2)
      n = 2
      call sub2(nums, n)
      write(6, *) 'nums after sub2:', nums(1), nums(2)
      call sub3(nums, n)
      write(6, *) 'nums after sub3:', nums(1), nums(2)
      end
