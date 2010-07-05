      program prog
      character s*127
      dimension numbers(20)
      do i=1,20
        numbers(i) = 253 + i
      enddo
      write(6, *) numbers(1), numbers(2), numbers(3), numbers(4),
     &  numbers(5), numbers(6), numbers(7), numbers(8)
      write(6, *) numbers(1), numbers(2), numbers(3), numbers(4),
     &  numbers(5), numbers(6), numbers(7), numbers(8), numbers(9),
     &  numbers(10)
      write(6, *) numbers(1), numbers(2), numbers(3), numbers(4),
     &  numbers(5), numbers(6), numbers(7), numbers(8), numbers(9),
     &  numbers(10), numbers(11), numbers(12)
      write(6, *) numbers(1), numbers(2), numbers(3), numbers(4),
     &  numbers(5), numbers(6), numbers(7), numbers(8), numbers(9),
     &  numbers(10), numbers(11), numbers(12), numbers(13), numbers(14)
      write(6, *) numbers(1), numbers(2), numbers(3), numbers(4),
     &  numbers(5), numbers(6), numbers(7), numbers(8), numbers(9),
     &  numbers(10), numbers(11), numbers(12), numbers(13), numbers(14),
     &  numbers(15), numbers(16)
      write(6, *) numbers(1), numbers(2), numbers(3), numbers(4),
     &  numbers(5), numbers(6), numbers(7), numbers(8), numbers(9),
     &  numbers(10), numbers(11), numbers(12), numbers(13), numbers(14),
     &  numbers(15), numbers(16), numbers(17), numbers(18)
      write(6, *) numbers(1), numbers(2), numbers(3), numbers(4),
     &  numbers(5), numbers(6), numbers(7), numbers(8), numbers(9),
     &  numbers(10), numbers(11), numbers(12), numbers(13), numbers(14),
     &  numbers(15), numbers(16), numbers(17), numbers(18), numbers(19),
     &  numbers(20)
      write(6, '(a)') '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklm
     &nopqrstuvwxyz)!@#$%^&*(`~-_+=[{]}\|;:''",<.>/?'
      write(6, '(a)') '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklm
     &nopqrstuvwxyz)!@#$%^&*\`~-_+=[{]}(|;:''",<.>/?'
      write(6, '(a)') '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklm
     &nopqrstuvwxyz)!@#$%^&*(\~-_+=[{]}`|;:''",<.>/?'
      s = 'qwertyuiopasdfghjklzxcvbnmqwertyuiopasdfghjklzxcvbnmqwertyuio
     &asdfghjkzxcvbnmqwerjkdfghjkertyjkxcghidfbndtyuiklmbvftyuiknbvdtyuh
     &'
      write(6, '(a)') s
 1    format(/,' Sorry, your unit cell, range of hkl, size of map,',
     &  ' and resolution will require ',/,
     &  ' redimensioning of the program.',//,
     &  ' This is quite easy:  You need to edit the source file ',
     &  ' for the program',/,
     &  ' and increase the value of "base_size" from ',i2,' to ',
     &  ' a larger value',//,
     &  '  Then recompile the program and try again.',
     &  //,' If you do not have the source code, then you can obtain',
     & / ' a version with a larger dimension from ',/,
     &  ' our web site.',/)
      write(6, 1) 12
 2    format(' first = ',F8.4,//,
     1  ' second:              ',F5.2,/,
     1  ' third:               ',F5.2,/,
     1  ' fourth:              ',F5.2,/,
     1  ' fifth:               ',F5.2,/,
     1  ' sixth:               ',F5.1)
      write(6, 2) 1.2, 3.4, 5.6, 7.8, 9.1, 2.3
      end
