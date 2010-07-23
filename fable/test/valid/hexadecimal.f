      program prog
      dimension nums(2)
      data nums /x'fe', x'dcba'/
      write(6, *) x'A'
      write(6, *) x'AB'
      write(6, *) x'ABC'
      write(6, *) x'ABCD'
      write(6, *) x'7FFFFFFF'
      write(6, *) nums
      end
