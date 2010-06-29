c1
      program prog!c2
c3
      dimension! c4
c5
     &nums(2)!c6
c7
      do i=1,2! c8
c9
        nums(i)!c10

c12
     &! c13
     &=!c14
     &i+47! c15
c16
      enddo!c17
      write(6, *, err=10) nums! c18
c19
      goto 20 !c20
c21
   10 stop 'write error' ! c22
c23
   20 continue!c24
c25
      end  !  c26
c27
