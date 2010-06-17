      program prog
      character*9 d
      character*8 t
      character*70 e
      call date(d)
      write(6, '(a)') d
      call time(t)
      write(6, '(a)') t
      call getenv(' PATH ', e)
      write(6, '(a)') e
      end
