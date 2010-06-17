      program prog
      character sa*2, sb*3, sc*5
      logical la, lb, lc
      lb = .not. la
      lc = la .and. lb
      lc = la .or. lb
      c = a + b
      c = a - b
      c = a * b
      c = a / b
      sa = 'x'
      sb = 'abc'
      sc = sa // sb
      end
