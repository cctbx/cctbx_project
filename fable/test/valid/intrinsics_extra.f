      program prog
      character*9 d
      character*8 t
      character*70 e
      character*8 c
      call date(d)
      write(6, '(a)') d
      call time(t)
      write(6, '(a)') t
      call getenv(' PATH ', e)
      write(6, '(a)') e
      call cpu_time(tm)
      write(6, '(f6.2)') tm
      c = 'echo YkD'
      i = system(c)
      write(6, *) i
      end
