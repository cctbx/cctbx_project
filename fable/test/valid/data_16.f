      parameter(num=2)
      real vals(num)
      data vals /num*1.2/
      write(6, *) (vals(i), i=1,num)
      end
