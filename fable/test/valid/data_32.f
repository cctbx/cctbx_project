      program prog
      dimension nums(2)
      dimension numsi(2)
      dimension numsj(2)
      character*2 str
      character*3 strs(2)
      character*4 strsi(2)
      character*5 strsj(2)
      data num /-34/
      data str /'YuIo'/
      data nums /+12, -34/
      data (numsi(i), i=1,2) /-56, +78/
      data strs /'Cde', 'FgH'/
      data (strsi(i), i=1,2) /'IjkL', 'MNOp'/
      data numj, numsj(1), strsj(1), numsj(2), strsj(2)
     &  /91, 23, 'Hjklo', 45, 'ASdfg'/
      write(6, '(i3)') num
      write(6, '(a)') str
      write(6, '(2i4)') nums
      write(6, '(2i4)') numsi
      write(6, '(2a)') strs
      write(6, '(2a)') strsi
      write(6, '(i2)') numj
      write(6, '(2i3)') numsj
      write(6, '(2a)') strsj
      call sub
      end

      subroutine sub
      dimension nums(2)
      data nums /-24,+35/
      write(6, '(2i4)') nums
      end
