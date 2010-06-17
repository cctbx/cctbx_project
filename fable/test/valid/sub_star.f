      subroutine sub1(num, *)
      if (num .eq. 1) return 1
      return
      end

      subroutine sub2(num, *, *)
      if (num .eq. 1) return 2
      if (num .eq. 2) return 1
      return
      end

      program prog
      do num=1, 2
        call sub1(num, *11)
        write(6, *) 'normal'
        goto 19
11      write(6, *) 'return 11'
19      continue
      enddo
      do num=1, 3
        call sub2(num, *21, *22)
        write(6, *) 'normal'
        goto 29
21      write(6, *) 'return 21'
        goto 29
22      write(6, *) 'return 22'
29      continue
      enddo
      end
