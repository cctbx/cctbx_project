      program prog
      character*1 first
      character*1 second
      character*2 third
      character*2 fourth
      first = 'x'
      second = 'y'
      first(1:1) = second(1:1)
      write(6, *) first
      third = 'pq'
      fourth = 'rs'
      third(1:2) = fourth(1:2)
      write(6, *) third
      third(2:2) = first(1:1)
      write(6, *) third
      third(1:1) = first
      write(6, *) third
      third(2:2) = "z"
      write(6, *) third
      call copy_overlapping
      call assign_directly
      end

      subroutine copy_overlapping
      character cval*26
      cval = 'abcdefghjiklmnopqrstuvwxyz'
      cval = cval(14:26)
      write(6, *) cval
      cval = 'abcdefghijklmnopqrtsuvwxyz'
      cval(5:11) = cval
      write(6, *) cval
      end

      subroutine assign_directly
      character small*2
      character big*3
      small = 'aB'
      big = small
      write(6, '(3a)') '[', big, ']'
      end
