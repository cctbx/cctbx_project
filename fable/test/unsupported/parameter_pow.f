      program prog
      save
      parameter(num1=2)
      parameter(num2=num1**2)
      dimension nums(num2)
      write(6, *) nums
      end
