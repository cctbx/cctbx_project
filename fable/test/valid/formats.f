      program prog

C "gerr" is not an error with ifort
C "ierr" was not tested with gfortran

      write(6, '(a)') 'ab'
      write(6, '(2a)') 'cd', 'ef'
      write(6, '(a3)') 'gh'
Cierr write(6, '(a3.1)') 'ij'
      write(6, '(2a3)') 'kl', 'mn'

Cgerr write(6, '(d)') 2.0
Cgerr write(6, '(2d)') 3.0, 4.0
Cierr write(6, '(d8)') 5.0
      write(6, '(d9.2)') 6.0
      write(6, '(2d9.2)') 7.0, 8.0

Cgerr write(6, '(e)') 2.0
Cgerr write(6, '(2e)') 3.0, 4.0
Cierr write(6, '(e8)') 5.0
      write(6, '(e9.2)') 6.0
      write(6, '(2e9.2)') 7.0, 8.0

Cgerr write(6, '(f)') 2.0
Cgerr write(6, '(2f)') 3.0, 4.0
Cierr write(6, '(f8)') 5.0
      write(6, '(f9.2)') 6.0
      write(6, '(2f9.2)') 7.0, 8.0

Cgerr write(6, '(g)') 2.0
Cgerr write(6, '(2g)') 3.0, 4.0
Cierr write(6, '(g8)') 5.0
      write(6, '(g9.2)') 6.0
      write(6, '(2g9.2)') 7.0, 8.0

Cierr write(6, '(ha)')
      write(6, '(1hb)')
      write(6, '(1hC)')
      write(6, '(2hdE)')
      write(6, '(3hF q)')

Cgerr write(6, '(i)') 2
Cgerr write(6, '(2i)') 3, 4
      write(6, '(i8)') 5
      write(6, '(i9.2)') 6
      write(6, '(2i9.2)') 7, 8
Cierr write(6, '(i+8)') 9

      write(6, '(l)') .false.
      write(6, '(2l)') .true., .false.
      write(6, '(l8)') .false.
Cierr write(6, '(l9.2)') .true.
      write(6, '(2l8)') .false., .true.

Cgerr write(6, '(z)') 17
Cgerr write(6, '(2z)') 18, 19
      write(6, '(z2)') 20
      write(6, '(z4.4)') 21
      write(6, '(z4.3)') 22
      write(6, '(2z4.3)') 23, 24

      write(6, '(a,t5,a)') 'ab', 'cd'
      write(6, '(a,tl2,a)') 'efg', 'hi'
      write(6, '(a,tr2,a)') 'jkl', 'mn'

      write(6, '(a,x,a)') 'abc', 'de'
      write(6, '(a,2x,a)') 'fgh', 'ij'

      write(6, '(a,/,a)') 'ab', 'cde'
      write(6, '(a/,a)') 'fg', 'hij'
      write(6, '(a,/a)') 'kl', 'mno'
      write(6, '(a/a)') 'pq', 'rst'

      write(6, '(a,:,a)') 'ab', 'cd'
      write(6, '(a,:,a)') 'ef'
      write(6, '(a,a)') 'gh'

      write(6, '(i2)') 2
      write(6, '(sp,i2)') 3
      write(6, '(spi2)') 4
      write(6, '(spi2,i2)') 5, 6
      write(6, '(sp,i2,ss,i2)') 7, 8
      write(6, '(sp,i2,ssi3)') 9, 10
      write(6, '(sp,i3,si3)') 11, 12
      write(6, '(sp,i3si3)') 13, 14 ! ifort ignores the s if comma is missing

      write(6, '(e9.2)') 2.0
Cgerr write(6, '(pe9.2)') 3.0
      write(6, '(1pe9.2)') 4.0
      write(6, '(1p,e9.2)') 5.0
      write(6, '(2pe9.2)') 6.0
      write(6, '(+2pe9.2)') 7.0
      write(6, '(-1pe9.2)') 8.0
      write(6, '(-2pe9.3)') 9.0
      write(6, '(-2p2e11.3)') 9.0, -2.0
      write(6, '(-2pe11.3,e12.4)') 10.0, -3.0

      write(6, '(bn,i3,bz,i2)') 2, 3
      write(6, '(bni3bzi2)') 4, 5

      write(6, '(i1$)') 2
      write(6, '(i1,$)') 3
      write(6, '(i1)') 4

      write(6, '(2(3hdef))')
      write(6, '(2(''def''))')

      end
