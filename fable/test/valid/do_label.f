      program prog
      do 10 num=1,2
        write(6, *) num
   10 continue
      do 20 num=1,2
        do 20 num2=3,4
          write(6, *) num, num2
   20 continue
      do 40 num=1,2
        do 30 num2=3,4
          write(6, *) num, num2
   30   write(6, *) num*10, num2*10
   40 continue
      end
