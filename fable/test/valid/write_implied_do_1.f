      program prog
      integer nums(2)
      nums(1) = 13
      nums(2) = -24
      write(6, '(2i4)') (nums(i), i=1,2)
      write(6, '(3i4)') (nums(i), (nums(j), j=1,2), i=1,2)
      write(6, '(5i4)') (
     &  nums(i),
     &    (nums(j), j=1,2),
     &    (nums(j), j=1,2),
     &      i=1,2)
      write(6, '(7i4)') (
     &  nums(i), (
     &    nums(j), (
     &      nums(k), k=1,2),
     &        j=1,2),
     &          i=1,2)
      end
