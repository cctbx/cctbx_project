      program prog
      parameter(num1=2)
      parameter(num2=num1+3)
      dimension nums(num2)
      do i=1,num2
        nums(i) = i*12
      enddo
      write(6, '(i2)') nums
      end
