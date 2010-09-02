      program examples
      implicit none
c
      byte b_var
c
      character s_var
      character*2 s2_var
c
      logical   l_var
      logical*1 l1_var
c
      integer   i_var
      integer*1 i1_var
      integer*2 i2_var
      integer*4 i4_var
      integer*8 i8_var
c
      real   r_var
      real*4 r4_var
      real*8 r8_var
      double precision d_var
      real*16 r16_var
c
      complex c_var
      complex*8 c8_var
      complex*16 c16_var
      double complex dc_var
      complex*32 c32_var
c
      open(unit=1, file='fable_tmp_d357d9d2', form='unformatted')
      write(1)
     & b_var,
     & s_var, s2_var,
     & l_var, l1_var,
     & i_var, i1_var, i2_var, i4_var, i8_var,
     & r_var, r4_var, r8_var, d_var, r16_var,
     & c_var, c8_var, c16_var, dc_var, c32_var
      close(1)
c
      write(6, '(a)') 'Done.'
c
      end
