      subroutine sub(nums_arg)
      parameter(isz=2)
      dimension nums_arg(isz)
      dimension nums_cmn(isz)
      common /scr/ nums_cmn
      write(6, *) nums_arg
      do i=1,2
        nums_cmn(i) = nums_cmn(i) + i * 19
      enddo
      write(6, *) nums_cmn
      end

      program prog
      dimension nums(2)
      do i=1,2
        nums(i) = 17*i
      enddo
      call sub(nums)
      call sub(nums)
      end
