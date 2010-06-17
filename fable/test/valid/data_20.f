      program prog
      dimension nums(2,3)
      data ((nums(i,j),j=1,3),i=1,2) /1,2,3,4,5,6/
      write(6, '(i1)') nums
      end
