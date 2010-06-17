      program prog
      dimension nums(2)
      character s2s(2)*2
      data (nums(i), s2s(i), i=1,2) /12, 'Xy', 34, 'Za'/
      do i=1,2
        write(6, '(i2, x, a)') nums(i), s2s(i)
      enddo
      end
