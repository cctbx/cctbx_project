      program prog
      character*(2) buf3
      character*(2) buf

      write(buf, '(a2)') "12"

      read(buf, fmt='(a3)') buf3
      write(*,*) '"', buf3, '"'

      end
