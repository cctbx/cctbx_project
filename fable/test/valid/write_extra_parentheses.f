      function i1(i)
      i1 = 7*i
      end

      function i2(i,j)
      i2 = i*8+j
      end

      function i3(i,j,k)
      i3 = i*9+j*29+k
      end

      program prog
      i = 3
      j = 4
      k = 5
      l = 6
      m = 7
      write(6, '(1i5)') -i*3
      write(6, '(1i5)') (-i*3)
      write(6, '(2i5)') (-i*3, -j*4)
      write(6, '(3i5)') (-i*3, (-j*4, -k*7))
      write(6, '(3i5)') ((-i*3, (-j*4, -k*7)))
      write(6, '(3i5)') ((-i*3, ((-j*4, -k*7))))
      write(6, '(4i5)') (-i*3, (-j*4, (-k*7, -l*8)))
      write(6, '(4i5)') ((-i*3, -j*4), (-k*7, -l*8))
      write(6, '(4i5)') (-i*3), -j*4, (-k*7, -l*8)
      write(6, '(5i5)') -i*3, (-j*4), -k*7, (-l*8, -m*9)
      write(6, '(4i5)') i1(-i*3), -j*4, (-k*7, -l*8)
      write(6, '(4i5)') i1(-i*3), (-j*4, (-k*7, -l*8))
      write(6, '(4i5)') i1(-i*3), ((-j*4, (-k*7, -l*8)))
      write(6, '(4i5)') (i1(-i*3)), ((-j*4, (-k*7, -l*8)))
      write(6, '(4i5)') ((i1(-i*3)), ((-j*4, (-k*7, -l*8))))
      write(6, '(2i5)') i1(-i*3), i3(-j*4, -k*7, -l*8)
      write(6, '(2i5)') (i1(-i*3), i3(-j*4, -k*7, -l*8))
      write(6, '(2i5)') ((i1(-i*3), i3(-j*4, -k*7, -l*8)))
      end
