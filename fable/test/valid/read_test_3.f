      program prog
      character*(2) buf3
      character*(3) buf33
      integer*4 k, l
      character*(2) buf

      write(buf, '(a2)') "12"

      read(buf, fmt='(a3)') buf3
      write(*,*) '"', buf3, '"'

      read(buf, fmt='(a3)') buf33
      write(*,*) '"', buf33, '"'

      read(buf, fmt='(i3,i3)', end=10) k, l
      write(*,*) k, l
      goto 20

 10   write(*,*) 'ERROR: in end'

 20   continue
      end
