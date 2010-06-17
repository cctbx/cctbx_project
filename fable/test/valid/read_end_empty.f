      program prog
      num = 0
  10  read(5, fmt='()', end=20)
      num = num + 1
      goto 10
  20  continue
      write(6, *) num
      end
