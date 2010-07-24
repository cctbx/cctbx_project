      program prog
      call exercise_write
      call exercise_read
      end

      subroutine exercise_write
      character c
      character s*2
      logical l
      integer i
      real r
      double precision d
      complex rc
      double complex dc
      open(1, file='fable_tmp_06c35b63',
     &  form='unformatted', status='unknown')
      c = 'x'
      write(1) c
      s = 'Ab'
      write(1) s
      l = .true.
      write(1) l
      i = 3
      write(1) i
      r = 4.1
      write(1) r
      d = 5.2D0
      write(1) d
      rc = cmplx(6.3, -7.4)
      write(1) rc
      dc = dcmplx(-8.5d0, 9.6d0)
      write(1) dc
      close(1)
      end

      subroutine exercise_read
      character c
      character s*2
      logical l
      integer i
      real r
      double precision d
      complex rc
      double complex dc
      open(1, file='fable_tmp_06c35b63',
     &  form='unformatted', status='old')
      read(1) c
      write(6, *) c
      read(1) s
      write(6, *) s
      read(1) l
      write(6, *) l
      read(1) i
      write(6, *) i
      read(1) r
      write(6, *) r
      read(1) d
      write(6, *) d
      read(1) rc
      write(6, *) rc
      read(1) dc
      write(6, *) dc
      close(1)
      end
