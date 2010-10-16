      program prog
      character s1*2
      character s2*3
      character s3*4
      s1 = 'x' // 'Y'
      write(6, '(a)') s1
      s1 = ('P') // 'q'
      write(6, '(a)') s1
      s1 = 'N' // ('o')
      write(6, '(a)') s1
      s3 = 'a' // 's' // 'd'
      write(6, '(a)') s3
      s3 = ('f' // 'g') // 'h'
      write(6, '(a)') s3
      s3 = 'j' // ('k' // 'l')
      write(6, '(a)') s3
      s2 = 'z' // s1
      write(6, '(a)') s2
      s2 = s1 // 'T'
      write(6, '(a)') s2
      s1 = 'U'
      s2 = s1 // 'v'
      write(6, '(a)') s2
      s2 = 'i'
      s3 = s1 // s2
      write(6, '(a)') s3
      s3 = s2 // s1
      write(6, '(a)') s3
      s1 = 'Wd'
      s3 = s1 // s2
      write(6, '(a)') s3
      s3 = s2 // s1
      write(6, '(a)') s3
      write(6, '(a)') 'Qwer' // 'ty'
      write(6, '(a)') ('Ui') // 'op'
      write(6, '(a)') 'M' // ('nb')
      write(6, '(a)') ('v' // ('cX')) // 'yz'
      end
