      program prog
      dimension nums(2)
      dimension data(2)
      data nums /12,34/
      data(1) = 1.2
      data(2) = 3.4
      write(6, '(f3.1,x,f3.1)') data
      write(6, '(i3,x,i3)', err=10) nums
   10 continue
      end
