      program prog
      character*(1) buf1
      character*(1) buf2
      character*(2) buf3
      character*(3) buf33
      character*(2) buf4(2)
      integer*4 k, l, m(2)
      double precision x
      integer ios, i
      character*(2) buf

      write(buf, '(a2)') "12"

      read(buf, fmt='(a1)', iostat=ios) buf1
      write(*,*) ios .lt. 0, ios .eq. 0, ios .gt. 0, '"', buf1, '"'

      read(buf, fmt='(a2)', iostat=ios) buf2
      write(*,*) ios .lt. 0, ios .eq. 0, ios .gt. 0, '"', buf2, '"'

      read(buf, fmt='(a3)', iostat=ios) buf3
      write(*,*) ios .lt. 0, ios .eq. 0, ios .gt. 0, '"', buf3, '"'

      read(buf, fmt='(a3)', iostat=ios) buf33
      write(*,*) ios .lt. 0, ios .eq. 0, ios .gt. 0, '"', buf33, '"'

      read(buf, fmt='(a3)', iostat=ios) (buf4(i), i=1,2)
      write(*,*) ios .lt. 0, ios .eq. 0, ios .gt. 0

      read(buf, fmt='(i3)', iostat=ios) k
      write(*,*) ios .lt. 0, ios .eq. 0, ios .gt. 0, k

      read(buf, fmt='(i3,i3)', iostat=ios) k, l
      write(*,*) ios .lt. 0, ios .eq. 0, ios .gt. 0, k, l

      read(buf, fmt='(i3)', iostat=ios) k, l
      write(*,*) ios .lt. 0, ios .eq. 0, ios .gt. 0

      read(buf, fmt='(i3)', iostat=ios) (m(i), i=1,2)
      write(*,*) ios .lt. 0, ios .eq. 0, ios .gt. 0

      read(buf, fmt='(d6.3)', iostat=ios) x
      write(*,*) ios .lt. 0, ios .eq. 0, ios .gt. 0, x

      write(buf, '(a2)') "0x"

      read(buf, fmt='(i2)', iostat=ios) k
      write(*,*) ios .lt. 0, ios .eq. 0, ios .gt. 0

      end
