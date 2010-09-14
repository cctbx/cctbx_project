      program prog
      character cvar*1024
      logical lvar
      open(10, file='fable_tmp_5d70aa2a')
      inquire(unit=10, name=cvar)
      i = istring_tail(cvar)
      if (i .lt. 18) then
        write(6, '(3a)') '[', 'FAILURE: ', cvar, ']'
      else
        write(6, '(3a)') '[', cvar(i-17:i), ']'
      endif
      inquire(10, exist=lvar)
      write(6, '(a,l1,a)') '(', lvar, ')'
      inquire(10, opened=lvar)
      write(6, '(a,l1,a)') '(', lvar, ')'
      inquire(file='fable_tmp_5d70aa2a', opened=lvar)
      write(6, '(a,l1,a)') '(', lvar, ')'
      inquire(unit=10, access=cvar)
      i = istring_tail(cvar)
      if (i .lt. 4) then
        write(6, '(3a)') '[', 'FAILURE: ', cvar, ']'
      else
        write(6, '(3a)') '[', cvar(i-3:i), ']'
      endif
      inquire(unit=10, sequential=cvar)
      i = istring_tail(cvar)
      if (i .lt. 2) then
        write(6, '(3a)') '[', 'FAILURE: ', cvar, ']'
      else
        write(6, '(3a)') '[', cvar(i-1:i), ']'
      endif
      inquire(unit=10, direct=cvar)
      i = istring_tail(cvar)
      if (i .lt. 2) then
        write(6, '(3a)') '[', 'FAILURE: ', cvar, ']'
      else
        write(6, '(3a)') '[', cvar(i-1:i), ']'
      endif
      inquire(unit=10, blank=cvar)
      i = istring_tail(cvar)
      if (i .lt. 4) then
        write(6, '(4a)') '[', 'FAILURE: ', cvar, ']'
      else
        write(6, '(3a)') '[', cvar(i-3:i), ']'
      endif
      write(10, '(a)') 'Abc'
      close(10)
      inquire(10, opened=lvar)
      write(6, '(a,l1,a)') '(', lvar, ')'
      inquire(file='fable_tmp_5d70aa2a', opened=lvar)
      write(6, '(a,l1,a)') '(', lvar, ')'
      inquire(file='fable_tmp_5d70aa2a', exist=lvar, err=10)
      write(6, '(a,l1,a)') '(', lvar, ')'
      goto 20
   10 write(6, '(a)') 'FAILURE inquire file'
   20 continue
      cvar = ' '
      inquire(file='fable_tmp_d185826b', name=cvar)
      i = istring_tail(cvar)
      if (i .lt. 18) then
        write(6, '(3a)') '[', 'FAILURE: ', cvar, ']'
      else
        write(6, '(3a)') '[', cvar(i-17:i), ']'
      endif
      end

      function istring_tail(s)
      character s*(*)
      do i=len(s),1,-1
        if (s(i:i) .ne. ' ') then
          istring_tail = i
          return
        endif
      enddo
      istring_tail = 0
      end
