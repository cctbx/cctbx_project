      program prog
      character*(3) buf33
      character*(2) buf4(2)
      integer*4 k, l, m(2)
      integer ios, i
      character*(2) buf

      write(buf, '(a2)') "12"

      read(buf, fmt='(a3)', iostat=ios, end=10, err=20) buf33
      write(*,*) ios .lt. 0, ios .eq. 0, ios .gt. 0, '"', buf33, '"'

      read(buf, fmt='(a3)', iostat=ios, end=101, err=20)
     @ (buf4(i), i=1,2)
      write(*,*) 'ERROR: 1'
      stop
 101  write(*,*) ios .lt. 0, ios .eq. 0, ios .gt. 0

      read(buf, fmt='(i3)', iostat=ios, end=102, err=20) k, l
      write(*,*) 'ERROR: 2'
      stop
 102  write(*,*) ios .lt. 0, ios .eq. 0, ios .gt. 0

      read(buf, fmt='(i3)', end=103, err=20, iostat=ios)
     @ (m(i), i=1,2)
      write(*,*) 'ERROR: 3'
      stop
 103  write(*,*) ios .lt. 0, ios .eq. 0, ios .gt. 0

      write(buf, '(a2)') "0x"

      read(buf, fmt='(i2)', iostat=ios, end=10, err=204) k
      write(*,*) 'ERROR: 4'
      stop
 204  write(*,*) ios .lt. 0, ios .eq. 0, ios .gt. 0

      write(*,*) 'OK'
      goto 30

 10   write(*,*) 'ERROR: in end='
      goto 30

 20   write(*,*) 'ERROR: in err='
      goto 30

 30   continue
      end
