      subroutine sub1(strs1)
      character*(*) strs1(*)
      do i=1,2
        write(6, '(a)') strs1(i)
      enddo
      end

      subroutine sub2(strs1)
      character*(*) strs1(*)
      do i=1,2
        strs1(i) = 'AbCd'
      enddo
      end

      program prog
      character strs2*4(2, 3)
      strs2(1,1) = 'A012'
      strs2(1,2) = 'C678'
      strs2(1,3) = 'E234'
      strs2(2,1) = 'B345'
      strs2(2,2) = 'D901'
      strs2(2,3) = 'F567'
      do j=1,3
        call sub1(strs2(1,j))
      enddo
      do j=1,3
        call sub2(strs2(1,j))
      enddo
      do j=1,3
        call sub1(strs2(1,j))
      enddo
      end
