
      logical function f(txt)
      character*(*) txt
      double precision value
      character*(1) buf

      open(10, file='fable_tmp_sio8miey', status='old')

      read(10, '(d1.0)', err=100, end=200) value
      write(*,*) 'value from file = ', value

      rewind(10)

      read(10, *, err=100, end=200) value
      write(*,*) 'value from file = ', value

      close(10)

      read(txt, *, err=100, end=200) value
      write(*,*) 'value from txt = ', value

      write(buf, '(a1,$)') '7'
      read(buf, '(d6.3)', err=100, end=200) value
      write(*,*) 'value from buf = ', value

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

      open(10, file='fable_tmp_sio8miey', status='unknown')
      write(10, '(a1,$)') '0'
      close(10)

      zero = '0'
      ff = f(zero)
      if(ff) write(*,*) 'OK'

      end
