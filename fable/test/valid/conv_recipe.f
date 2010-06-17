C cctbx_project/compcomm/newsletter09/conv_recipe.py, svn rev. 9983

      subroutine show_resolution(h, k, l)
      implicit none
      integer h, k, l
      real a, b, c
      logical first
      real ass, bss, css
      real dss
      COMMON /abc/ a, b, c
      SAVE first
      SAVE ass, bss, css
      DATA first /.true./
      if (first) then
        first = .false.
        if (a .le. 0 .or. b .le. 0 .or. c .le. 0) then
          write(0, '(1x,a)') 'invalid unit cell constants.'
          stop
        endif
        ass = 1/(a*a)
        bss = 1/(b*b)
        css = 1/(c*c)
      endif
      dss = h*h*ass + k*k*bss + l*l*css
      if (dss .eq. 0) then
        write(6, '(3(1x,i3),1x,a)')
     &    h, k, l, '    infinity'
      else
        write(6, '(3(1x,i3),1x,f12.6)')
     &    h, k, l, sqrt(1/dss)
      endif
      return
      end

      PROGRAM conv_recipe
      implicit none
      real a, b, c
      COMMON /abc/ a, b, c
      a = 11.0
      b = 12.0
      c = 13.0
      call show_resolution(0, 0, 0)
      call show_resolution(1, 2, 3)
      end
