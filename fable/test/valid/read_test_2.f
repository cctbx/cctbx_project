      program prog
      integer*4 k, l
      character*(2) buf

      write(buf, '(a2)') "12"

      read(buf, fmt='(i3,i4)', end=10) k, l
      write(*,*) k, l
      write(*,*) 'OK'
      goto 20

 10   write(*,*) 'ERROR: in end'

 20   continue
      end
