      program prog
      dimension nums(4)
      data (nums(i), i=1,4,2) /12,34/
      data (nums(i), i=2,4,2) /56,78/
      write(6, '(4i3)') nums
      end
