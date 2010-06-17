      subroutine sub(sz, nums)
      save
      integer sz
      dimension nums(sz)
      do i=1,sz
        write(6, *) nums(i)
      enddo
      end

      program prog
      dimension nums(2)
      call sub(2, nums)
      end
