      program prog
      external ext_impl
      call sub(ext_impl)
      end

      subroutine sub(ext)
      external ext
      call sub2(ext)
      call ext(3)
      end

      subroutine sub2(ext)
      external ext
      call ext(4)
      end

      subroutine ext_impl(num)
      write(6, *) num
      end
