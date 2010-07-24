      program prog
      character sa*2, sb*3
      logical la, lb
      lb = .not. la
      write(6, *) la .and. lb
      write(6, *) la .or. lb
      b = 1
      write(6, *) a + b
      write(6, *) a - b
      write(6, *) a * b
      write(6, *) a / b
      sa = 'x'
      sb = 'abc'
      write(6, *) sa // sb
      write(6, *) a .eq. b
      write(6, *) a .ne. b
      write(6, *) a .lt. b
      write(6, *) a .le. b
      write(6, *) a .gt. b
      write(6, *) a .ge. b
      write(6, *) a == b
      write(6, *) a /= b
      write(6, *) a < b
      write(6, *) a <= b
      write(6, *) a > b
      write(6, *) a >= b
      end
