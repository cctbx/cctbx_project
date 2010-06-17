      subroutine sub(str)
      character str*(*)
      if (str(1:1) .eq. ' ') then
        write(6, '(a)') 'str starts with a blank'
      else
        write(6, '(a)') 'str does not start with a blank'
      endif
      end

      program prog
      character str2*2
      str2 = ' '
      call sub(str2)
      str2 = 'x'
      call sub(str2)
      end
