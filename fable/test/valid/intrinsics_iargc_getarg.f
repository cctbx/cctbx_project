      subroutine sub1(n, buf)
      character buf*(*)(*)
      n = iargc()
      do i=1,n
        call getarg(i, buf(i))
      enddo
      end

      subroutine sub2(n, buf)
      implicit none
      integer n
      character buf*(*)(*)
      integer i, j
      n = iargc()
      j = 1
      do i=n,1,-1
        call getarg(i, buf(j))
        j = j + 1
      enddo
      end

      subroutine sub3(n, buf)
      implicit none
      integer n
      character buf*(*)(*)
      integer i
      integer iargc
      n = iargc()
      do i=2,n
        call getarg(i, buf(i-1))
      enddo
      if (iargc() .ne. 0) then
        call getarg(1, buf(n))
      endif
      end

      program prog
      character*4 buf(10)
      call sub1(n, buf)
      write(6, '(a)') 'A', (buf(i), i=1,n)
      call sub2(n, buf)
      write(6, '(a)') 'B', (buf(i), i=1,n)
      call sub3(n, buf)
      write(6, '(a)') 'C', (buf(i), i=1,n)
      end
