      program prog
      common /first/ n1, n2(2)
      dimension nums(3)
      equivalence(nums, n1)
      equivalence(m2, n2)
      n1 = 25
      n2(1) = 93
      n2(2) = 37
      write(6, *) nums
      call sub
      write(6, *) nums
      write(6, *) m2
      end

      subroutine sub
      common /second/ n2
      n2 = 54
      write(6, *) n2
      end
