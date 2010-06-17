      subroutine sub1(str)
      character str*(*)
      write(6, '(i1)') len(str)
      write(6, '(a)') str
      end

      subroutine sub2(small)
      character small*(*)
      character big*2
      big = small
      write(6, '(3a)') '[', big, ']'
      end

      program prog
      character str2*2
      character str3*3
      str2 = 'Pq'
      str3 = 'rSt'
      call sub1(str2)
      call sub1(str3)
      call sub2('a')
      end
