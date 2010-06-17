      subroutine sub1
      dimension nums(1)
      equivalence(nums(1), n1)
      nums(1) = 12
      write(6, *) n1
      end

      subroutine sub2
      dimension nums(2)
      equivalence(nums(1), n1)
      equivalence(nums(2), n2)
      nums(1) = 34
      nums(2) = 56
      write(6, *) n1, n2
      end

      subroutine sub3
      dimension nums(3)
      equivalence(nums(1), n1)
      equivalence(nums(2), n2)
      equivalence(nums(3), n3)
      nums(1) = 78
      nums(2) = 90
      nums(3) = 23
      write(6, *) n1, n2, n3
      end

      subroutine sub4
      dimension nums(4)
      equivalence(nums(1), n1)
      equivalence(nums(2), n2)
      equivalence(nums(3), n3)
      equivalence(nums(4), n4)
      nums(1) = 45
      nums(2) = 67
      nums(3) = 89
      nums(4) = 43
      write(6, *) n1, n2, n3, n4
      end

      subroutine sub5
      dimension nums(5)
      equivalence(nums(1), n1)
      equivalence(nums(2), n2)
      equivalence(nums(3), n3)
      equivalence(nums(4), n4)
      equivalence(nums(5), n5)
      nums(1) = 65
      nums(2) = 87
      nums(3) = 89
      nums(4) = 91
      nums(5) = 93
      write(6, *) n1, n2, n3, n4, n5
      end

      subroutine sub6
      dimension nums(6)
      equivalence(nums(1), n1)
      equivalence(nums(2), n2)
      equivalence(nums(3), n3)
      equivalence(nums(4), n4)
      equivalence(nums(5), n5)
      equivalence(nums(6), n6)
      nums(1) = 94
      nums(2) = 95
      nums(3) = 96
      nums(4) = 97
      nums(5) = 27
      nums(6) = 29
      write(6, *) n1, n2, n3, n4, n5, n6
      end

      program prog
      call sub1
      call sub2
      call sub3
      call sub4
      call sub5
      call sub6
      end
