      integer function fun(num)
      write(6, *) num
      fun = 10 + num
      return
      end

      subroutine sub1(func)
      implicit none
      external func
      integer func
      integer i
      save
      i = func(i)
      write(6, *) i
      return
      end

      subroutine sub2(func)
      external func
      integer func
      i = func(3)
      write(6, *) i
      return
      end

      program prog
      implicit none
      external fun
      integer fun
      integer i
      i = fun(1)
      write(6, *) i
      call sub1(fun)
      call sub1(fun)
      call sub2(fun)
      end
