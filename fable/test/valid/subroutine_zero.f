      subroutine zero(vals)
      dimension vals(2)
      vals(1) = 12
      vals(2) = 34
      end

      program prog
      dimension vals(2)
      call zero(vals)
      write(6, '(2f5.1)') vals(1), vals(2)
      end
