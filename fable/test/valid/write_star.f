      program prog
      character c
      character s1*1, s2*2
      logical l
      integer i
      integer*8 j
      real r
      double precision d
      c = 'x'
      write(6,*) c
      write(6,*) 'i is zero.'
      write(6,*) l
      write(6,*) i
      write(6,*) j
      write(6,*) r
      write(6,*) d
      write(6,*), 1.d111
      write(6,*) -1.d111
      write(6,*)
      write(6,*) c, c, c, c, c, c, c, c, c, c, c, c
      write(6,*) 'i is ', 'zero', '.'
      write(6,*) l, l, l, l, l, l, l, l, l, l, l, l
      write(6,*) i, i, i, i, i, i, i, i, i, i, i, i
      write(6,*) j, j, j, j, j, j, j, j, j, j, j, j
      write(6,*) r, r, r, r, r, r, r, r, r, r, r, r
      write(6,*) d, d, d, d, d, d, d, d, d, d, d, d
      s1 = 'x'
      s2 = 'yz'
      write(6,*) s1, s1
      write(6,*) s1, s2
      write(6,*) s2, s1
      write(6,*) s2, s2
      write(6,*) s1, 12
      write(6,*) s2, 34
      write(6,*) 56, s1
      write(6,*) 78, s2
      write(6,*) 'aBcD ', 12, ' eFgHi ', 345
      end
