      program test
      character*(3) buf33
      character*(1) buf
      integer*4 ios
      write(buf, '(a1,$)') '7'
      read(buf, '(a3)', err=100, end=200, iostat=ios) buf33
      write(*,*) 'buf33 from buf = ', '"',buf33,'"', ' ios =', ios
      goto 999
 100  print *, 'ERROR: err'
      goto 999
 200  print *, 'ERROR: end'
      goto 999
 999  continue
      end
