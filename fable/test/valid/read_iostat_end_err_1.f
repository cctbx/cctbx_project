
      logical function f(txt)
      character*(*) txt
      double precision value
      character*(3) buf33
      character*(1) buf
      integer*4 ios

      open(10, file='fable_tmp_sa1oosh4', status='old')

      read(10, '(d1.0)', err=100, end=200, iostat=ios) value
      write(*,*) 'value from file = ', value, 'ios =', ios

      rewind(10)

      read(10, *, err=100, end=200, iostat=ios) value
      write(*,*) 'value from file = ', value, 'ios =', ios

      close(10)

      read(txt, *, err=100, end=200, iostat=ios) value
      write(*,*) 'value from txt = ', value, 'ios =', ios

      write(buf, '(a1,$)') '7'

      read(buf, '(d6.3)', err=100, end=200, iostat=ios) value
      write(*,*) 'value from buf = ', value, 'ios =', ios

      read(buf, '(a3)', err=100, end=200, iostat=ios) buf33
      write(*,*) 'buf33 from buf = ', '"',buf33,'"', ' ios =', ios


      f = .true.
      return

 100  continue
      print *, 'ERROR: err'
      f = .false.
      return

 200  continue
      print *, 'ERROR: end'
      f = .false.
      return

      end

      program read_err_end
      logical f, ff
      character*(1) zero

      open(10, file='fable_tmp_sa1oosh4', status='unknown')
      write(10, '(a1,$)') '0'
      close(10)

      zero = '0'
      ff = f(zero)
      if(ff) write(*,*) 'OK'

      end
