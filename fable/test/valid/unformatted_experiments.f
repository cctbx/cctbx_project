      subroutine try_integers
      open(unit=1, form='unformatted', status='unknown')
      write(1) 12, 34
      write(1) 56
      write(1) 78, 90
      close(1)
      open(unit=1, form='unformatted', status='old')
      read(1) num1, num2
      write(6, *) num1, num2
      read(1) num1
      write(6, *) num1
      read(1) num1, num2
      write(6, *) num1, num2
      close(1)
      open(unit=1, form='unformatted', status='old')
      read(1) num1
      write(6, *) num1
      read(1) num1
      write(6, *) num1
      read(1) num1
      write(6, *) num1
      close(1)
      end

      subroutine try_strings
      character str2*2, str3*3, str8*8
      open(unit=1, form='unformatted', status='unknown')
      write(1) 'Ab', 'cDe'
      write(1) 'f', 'Ghi'
      write(1) 'JkLmnop', 'qrstuvw'
      close(1)
      open(unit=1, form='unformatted', status='old')
      read(1) str2, str3
      write(6, *) '[', str2, '], [', str3, ']'
      read(1) str2
      write(6, *) '[', str2, ']'
      read(1) str2, str3
      write(6, *) '[', str2, '], [', str3, ']'
      close(1)
      open(unit=1, form='unformatted', status='old')
      read(1) str2
      write(6, *) '[', str2, ']'
      read(1) str2
      write(6, *) '[', str2, ']'
      read(1) str2
      write(6, *) '[', str2, ']'
      close(1)
      if (.false.) then
        open(unit=1, form='unformatted', status='old')
        read(1) str8
C         ifort: input statement requires too much data
C         gfortran: Short record on unformatted read
        write(6, *) '[', str8, ']'
        close(1)
      endif
      end

      subroutine try_mixed_integers
      integer*4 i4o, i4a, i4b
      integer*8 i8o, i8a, i8b
      open(unit=1, form='unformatted', status='unknown')
      i4o = 4321
      i8o = 1234
      i8o = i8o * 2**16
      i8o = i8o * 2**16
      i8o = i8o + 5678
      write(1) i4o, i4o+1111
      close(1)
      open(unit=1, form='unformatted', status='old')
      read(1) i4a, i4b
      write(6, *) i4a, i4b
      close(1)
      open(unit=1, form='unformatted', status='old')
      read(1) i8a
      write(6, *) i8a
      close(1)
      open(unit=1, form='unformatted', status='unknown')
      write(1) i8o
      close(1)
      open(unit=1, form='unformatted', status='old')
      read(1) i8a
      write(6, *) i8a
      close(1)
      open(unit=1, form='unformatted', status='old')
      read(1) i4a, i4b
      write(6, *) i4a, i4b
      close(1)
      end

      program prog
      call try_integers
      call try_strings
      call try_mixed_integers
      end
