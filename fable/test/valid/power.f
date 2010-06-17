      program prog
      parameter(val=1.2**3.4)
      write(6, '(f5.3)') val
      x = 1.2 + 3.4**5.6 / 7.8
      write(6, '(f5.1)') x
      x = (1.2 + 3.4)**5.6 / 7.8
      write(6, '(f5.1)') x
      x = (1.2 + 3.4)**(5.6 / 7.8)
      write(6, '(f5.3)') x
      x = -1.3**2
      write(6, '(f5.2)') x
      x = (-1.3)**2
      write(6, '(f4.2)') x
      end
