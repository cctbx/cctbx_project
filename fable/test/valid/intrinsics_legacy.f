      program intrinsics_2
      integer*4 k,l,m,s,i,j
      character*(2) ss(6)
      double precision xx
      real*4 x
      k = -1
      l = 0
      m = 1
      s = -1
      write(*,*) isign(k,s), isign(l,s), isign(m,s)
      s = 0
      write(*,*) isign(k,s), isign(l,s), isign(m,s)
      s = 1
      write(*,*) isign(k,s), isign(l,s), isign(m,s)
      xx = 7.2973525698D-3
      x = sngl(xx)
      write(*,*) xx, x, dble(x)
      write(*,*) idnint(-1.37D0), idnint(-1.73D0)
      write(*,*) idnint(xx), idnint(0.5D0), idnint(1.5D0)
      write(*,*) idnint(1.37D0), idnint(1.73D0)
      ss(1) = 'ab'
      ss(2) = 'AB'
      ss(3) = '0x'
      ss(4) = ',!'
      ss(5) = ' b'
      ss(5) = '# '
      write(*,*) ((llt(ss(i),ss(j)) , i=1,6), j=1,6)
      write(*,*) ((lle(ss(i),ss(j)) , i=1,6), j=1,6)
      write(*,*) ((lgt(ss(i),ss(j)) , i=1,6), j=1,6)
      write(*,*) ((lge(ss(i),ss(j)) , i=1,6), j=1,6)
      end
      
