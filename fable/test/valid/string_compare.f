      program prog
      character str1*1, str2*2
      str1 = ' '
      str2 = ' '
      write(6, '(a,2l1)') 'a', (str1 .eq. str2), (str1 .ne. str2)
      str1 = 'x'
      write(6, '(a,2l1)') 'b', (str1 .eq. str2), (str1 .ne. str2)
      str2 = 'x'
      write(6, '(a,2l1)') 'c', (str1 .eq. str2), (str1 .ne. str2)
      str2 = 'xy'
      write(6, '(a,2l1)') 'd', (str1 .eq. str2), (str1 .ne. str2)
      str2 = ' y'
      write(6, '(a,2l1)') 'e', (str1 .eq. str2), (str1 .ne. str2)
      str1 = ' '
      write(6, '(a,2l1)') 'f', (str1 .eq. str2), (str1 .ne. str2)

      str1 = ' '
      write(6, '(a,2l1)') 'g', (str1 .eq. ' '), (str1 .ne. ' ')
      str1 = 'x'
      write(6, '(a,2l1)') 'h', (str1 .eq. 'x'), (str1 .ne. 'x')
      write(6, '(a,2l1)') 'i', (str1 .eq. 'x '), (str1 .ne. 'x ')
      write(6, '(a,2l1)') 'j', (str1 .eq. 'xy'), (str1 .ne. 'xy')

      str2 = ' '
      write(6, '(a,2l1)') 'k', (str2 .eq. ' '), (str2 .ne. ' ')
      str2 = 'x'
      write(6, '(a,2l1)') 'l', (str2 .eq. 'x'), (str2 .ne. 'x')
      write(6, '(a,2l1)') 'm', (str2 .eq. 'x '), (str2 .ne. 'x ')
      write(6, '(a,2l1)') 'n', (str2 .eq. 'xy'), (str2 .ne. 'xy')
      str2 = 'xy'
      write(6, '(a,2l1)') 'o', (str2 .eq. 'xy'), (str2 .ne. 'xy')
      str2 = ' y'
      write(6, '(a,2l1)') 'p', (str2 .eq. ' y'), (str2 .ne. ' y')
      write(6, '(a,2l1)') 'q', (str2 .eq. ' y '), (str2 .ne. ' y ')
      write(6, '(a,2l1)') 'r', (str2 .eq. ' yz'), (str2 .ne. ' yz')
      write(6, '(a,2l1)') 's', (' y' .eq. str2), (' y' .ne. str2)

      call lt_le_gt_ge
      end

      subroutine lt_le_gt_ge
      character ab*2, abc*3, ab_*3, ac*2
      ab = 'ab'
      abc = 'abc'
      ab_ = 'ab '
      ac = 'ac'
      write(6, '(a)') 'ab ab'
      write(6, *) ab .lt. ab
      write(6, *) ab .le. ab
      write(6, *) ab .gt. ab
      write(6, *) ab .ge. ab
      write(6, '(a)') 'ab abc'
      write(6, *) ab .lt. abc
      write(6, *) ab .le. abc
      write(6, *) ab .gt. abc
      write(6, *) ab .ge. abc
      write(6, '(a)') 'abc ab'
      write(6, *) abc .lt. ab
      write(6, *) abc .le. ab
      write(6, *) abc .gt. ab
      write(6, *) abc .ge. ab
      write(6, '(a)') 'ab ab_'
      write(6, *) ab .lt. ab_
      write(6, *) ab .le. ab_
      write(6, *) ab .gt. ab_
      write(6, *) ab .ge. ab_
      write(6, '(a)') 'ab_ ab'
      write(6, *) ab_ .lt. ab
      write(6, *) ab_ .le. ab
      write(6, *) ab_ .gt. ab
      write(6, *) ab_ .ge. ab
      write(6, '(a)') 'ab ac'
      write(6, *) ab .lt. ac
      write(6, *) ab .le. ac
      write(6, *) ab .gt. ac
      write(6, *) ab .ge. ac
      write(6, '(a)') 'ac ab'
      write(6, *) ac .lt. ab
      write(6, *) ac .le. ab
      write(6, *) ac .gt. ab
      write(6, *) ac .ge. ab
      end
