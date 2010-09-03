      program prog
      call exercise_misc
      call exercise_wrap
      call exercise_logical
      call exercise_complex
      call exercise_escape_string_literal
      end

      subroutine exercise_misc
      num = 12
      val = 23.45
      write(6, '(i3)') num
      write(6, '(f6.2)') val
      write(6, '(i3,1x,f6.2)') num, val
      write(6, '(2i3)') num, num+1
      write(6, '(1(i3))') num
      write(6, '(2(i3))') num, num+1
      write(6, '(1(1(i3)))') num
      write(6, '(1(2(i3)))') num, num+1
      write(6, '(2(1(i3)))') num, num+1
      write(6, '(2(2i3))') num, num+1, num+2, num+3
      write(6, '(2(i3,1x,f6.2))') num, val, num+1, val+1
      write(6, '(2(i3,1x,2(f6.2,2x,i2),''_''))')
     &  num, val, num+1, val+1, num+2,
     &  num+3, val+2, num+4, val+3, num+5
      write(6, '(i2)') 1,2
      write(6, '(i2,i3)') 1,2,3
      write(6, '(i2,i3,''_'')') 1,2,3
      write(6, '(i2,i3,''_'')') 1,2,3,4
      write(6, '(2(i2,i3,''_''))') 1,2,3,4,5
      write(6, '(2(i2,i3,''_''))') 1,2,3,4,5,6
      write(6, '(''>'',2(i2,i3,''_''))') 1,2,3,4,5,6
      write(6, '(i4,2(i2,i3,''_''))') -1,1,2,3,4,5,6
      write(6, '(i2,2(i3,2(i4,i5)))')
     &  1, 2,3,4,5,6, 7,8,9,10,11, 12
      write(6, '(i2,2(i3,2(i4,i5),i6))')
     &  1, 2,3,4,5,6,7, 8,9,10,11,12,13, 14
      write(6, '(i2,2(i2,i3),2(i4,i5))') 1,2,3,4,5,6,7,8,9,10
      write(6, '(i2$)') 64
      write(6, '(i2)') 53
      write(6, '(2(i2,(3i2)))') 1,2,3,4,5,6,7,8
      write(6, '((2i2,(3i2)))') 1,2,3,4,5
c
      write(6, '(2(i2,1x))') 12, 34
      write(6, '(2(i2,'' ''))') 12, 34
      write(6, '(2(i2,''A''))') 12, 34
      end

      subroutine exercise_wrap
      write(6, '(i2)') 1,2,3,4
      write(6, '((i2))') 1,2,3,4
      write(6, '(((i2)))') 1,2,3,4
      write(6, '(2i2)') 1,2,3,4
      write(6, '((2i2))') 1,2,3,4
      write(6, '(((2i2)))') 1,2,3,4
      write(6, '(i2,i3)') 1,2,3,4
      write(6, '(i2,i3)') 1,2,3,4,5
      write(6, '((i2,i3))') 1,2,3,4
      write(6, '((i2,i3))') 1,2,3,4,5
      write(6, '(((i2,i3)))') 1,2,3,4
      write(6, '(((i2,i3)))') 1,2,3,4,5
      write(6, '(2i2,i3)') 1,2,3,4,5,6
      write(6, '(2i2,i3)') 1,2,3,4,5,6,7
      write(6, '(i2,2i3)') 1,2,3,4,5,6
      write(6, '(i2,2i3)') 1,2,3,4,5,6,7
      write(6, '((i2),i3)') 1,2,3,4
      write(6, '(i2,(i3))') 1,2,3,4,5
      write(6, '((i2,i3),i4)') 1,2,3,4,5,6
      write(6, '((i2,i3),i4)') 1,2,3,4,5,6,7
      write(6, '((i2,i3),i4)') 1,2,3,4,5,6,7,8
      write(6, '(i2,(i3,i4))') 1,2,3,4,5,6
      write(6, '(i2,(i3,i4))') 1,2,3,4,5,6,7
      write(6, '(i2,(i3,i4))') 1,2,3,4,5,6,7,8
      write(6, '((i2,(i3,(i4))))') 1,2,3,4,5,6
      write(6, '((i2,(i3,(i4))))') 1,2,3,4,5,6,7
      write(6, '((i2,(i3,(i4))))') 1,2,3,4,5,6,7,8
      write(6, '((i2,(i3,(i4))))') 1,2,3,4,5,6,7,8,9
      write(6, '((i2,(i3,(i4))))') 1,2,3,4,5,6,7,8,9,10
      write(6, '(2(i2,i3))') 11,12,21,22,23,24
      write(6, '(i2,2(i3,i4),i5)') 1,2,3,4,5,6,7,8,9,10,11,12,13
      write(6, '(i2,2(i3,i4),(i5,i6))') 1,2,3,4,5,6,7,8,9,10,11,12
      write(6, '(i2,2(i3,i4),2(i5,i6))') 1,2,3,4,5,6,7,8,9,10,11,12
      write(6, '(i2,(i3,i4),2(i5,i6))') 1,2,3,4,5,6,7,8,9,10,11,12
      write(6, '((2i2,(3i3)))') 1,2,3,4,5,6,7,8,9,10,11,12,13,14
      end

      subroutine exercise_logical
      logical flag
      do i=1,2
        flag = i .eq. 2
        write(6, '(l1)') flag
        write(6, '(l2)') flag
        write(6, '(l3)') flag
        write(6, '(l4)') flag
        write(6, '(l5)') flag
        write(6, '(l6)') flag
      enddo
      write(6, '(a,l1)') 'Abc', .false.
      write(6, '(a,1x,l1)') 'Xyz', .true.
      do i=1,2
        flag = i .eq. 2
        write(6, '(a,5x,l2)') 'Ab', flag
        write(6, '(a,4x,l3)') 'cD', flag
        write(6, '(a,3x,l4)') 'EF', flag
        write(6, '(a,2x,l5)') 'gh', flag
        write(6, '(a,1x,l6)') 'iJ', flag
      enddo
      write(6, '(a,l1)') 'Abc', .false., 'pQr'
      write(6, '(a,1x,l1,a)') 'Xyz', .true., 'tuV'
      do i=1,2
        flag = i .eq. 2
        write(6, '(a,5x,l2,a)') 'Ab', flag, 'kl'
        write(6, '(a,4x,l3,a)') 'cD', flag, 'mN'
        write(6, '(a,3x,l4,a)') 'EF', flag, 'OP'
        write(6, '(a,2x,l5,a)') 'gh', flag, 'Qp'
        write(6, '(a,1x,l6,a)') 'iJ', flag, 'rS'
      enddo
      end

      subroutine exercise_complex
      parameter(n=4)
      dimension w(n), a(2,n), b(2,n), c(n), d(n), e(n)
      do i=1,n
        w(i) = 10+i
        do j=1,2
          a(j,i) = 20+i+j/10.
          b(j,i) = 30+i+j/10.
        enddo
        c(i) = 40+i
        d(i) = 50+i
        e(i) = 60+i
      enddo
      write(6, 10)
     &  n,
     &  (i, w(i),
     &    (a(k1,i),k1=1,2),
     &    (b(kk1,i),kk1=1,2),
     &    c(i), d(i), e(i),
     &  i=1, n)
 10   format(/, 'Ab', i2, 'Cd', 'Ef', /, 'Gh', /,
     &  3(/, i3, 'Ij', f8.4, /, 'Kl', 'Mn', /,
     &    ('Op', 4f7.3, 3f7.3, /)))
      end

      subroutine exercise_escape_string_literal
      write(6, '(a)') '\"	??'
      end
