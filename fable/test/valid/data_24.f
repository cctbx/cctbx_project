      program prog
      dimension nums(2,2)
      character s2s(2,2)*2
      data ((nums(i,j), s2s(i,j), i=1,2,2), j=1,2),
     &     ((nums(i,j), s2s(i,j), i=2,2,2), j=1,2)
     &  /12, 'Xy', 34, 'Za', 56, 'cD', 78, 'eF'/
      do i=1,2
        do j=1,2
          write(6, '(i2, x, a)') nums(i,j), s2s(i,j)
        enddo
      enddo
      end
