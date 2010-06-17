      program prog
      num = 84
      val = float(num)
      val = val + 0.99
      num = int(val)
      write(6, '(f4.1,1x,i2)') val, num
      end
