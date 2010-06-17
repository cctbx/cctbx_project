      program prog
      character*3 s3(2)
      character*8 s8
      common /com/ s3, s8
      character*6 s6
      equivalence(s3, s6)
      character*2 s2(3)
      equivalence(s3, s2)
      character*4 s4(2)
      equivalence(s8, s4)
      character*8 s8e
      equivalence(s8, s8e)
      character*1 s1(5)
      equivalence(s1, s3(2))
      s3(1) = 'AbC'
      s3(2) = 'dEf'
      write(6, '(a)') s6
      s6 = 'PqrStu'
      do i=1,3
        write(6, '(a)') s2(i)
      enddo
      s2(1) = 'gH'
      s2(2) = 'Ij'
      s2(3) = 'kL'
      do i=1,2
        write(6, '(a)') s3(i)
      enddo
      s8e = 'AbcdefgH'
      do i=1,2
        write(6, '(a)') s4(i)
      enddo
      do i=1,5
        write(6, '(a)') s1(i)
      enddo
      end
