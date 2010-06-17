C cctbx_project/compcomm/newsletter09/sf_times.py, svn rev. 9983
C cos_wrapper, exp_wrapper replaced with COS, SIN, EXP

      subroutine sf(abcss, n_scatt, xyz, b_iso, n_refl, hkl, f_calc)
      implicit none
      REAL abcss(3)
      integer n_scatt
      REAL xyz(3, *)
      REAL b_iso(*)
      integer n_refl
      integer hkl(3, *)
      REAL f_calc(2, *)
      integer i_refl, i_scatt, j, h
      REAL phi, cphi, sphi, dss, ldw, dw, a, b
      DO i_refl=1,n_refl
        a = 0
        b = 0
        DO i_scatt=1,n_scatt
          phi = 0
          DO j=1,3
            phi = phi + hkl(j,i_refl) * xyz(j,i_scatt)
          enddo
          phi = phi * 2 * 3.1415926535897931
          cphi = COS(phi)
          sphi = SIN(phi)
          dss = 0
          DO j=1,3
            h = hkl(j,i_refl)
            dss = dss + h*h * abcss(j)
          enddo
          ldw = -0.25 * dss * b_iso(i_scatt)
          dw = EXP(ldw)
          a = a + dw * cphi
          b = b + dw * sphi
        enddo
        f_calc(1, i_refl) = a
        f_calc(2, i_refl) = b
      enddo
      return
      end

      program run
      implicit none
      REAL abcss(3)
      integer n_scatt
      parameter(n_scatt=10)
      REAL xyz(3, n_scatt)
      REAL b_iso(n_scatt)
      integer n_refl
      parameter(n_refl=10)
      integer hkl(3, n_refl)
      REAL f_calc(2, n_refl)
      integer i, j, jr
      REAL a, b, max_a, max_b
      abcss(1) = 1/(11.0*11.0)
      abcss(2) = 1/(12.0*12.0)
      abcss(3) = 1/(13.0*13.0)
      jr = 0
      DO i=1,n_scatt
        DO j=1,3
          jr = mod(jr*1366+150889, 714025)
          xyz(j,i) = (mod(jr, 20000) - 10000) / 10000.0
        enddo
      enddo
      DO i=1,n_scatt
        jr = mod(jr*1366+150889, 714025)
        b_iso(i) = mod(jr, 10000) / 100.0
      enddo
      if (n_scatt .le. 10) then
        DO i=1,n_scatt
          write(6, '(4(1x,f9.6))')
     &      xyz(1,i), xyz(2,i), xyz(3, i), b_iso(i)
        enddo
      endif
      DO i=1,n_refl
        DO j=1,3
          jr = mod(jr*1366+150889, 714025)
          hkl(j,i) = mod(jr, 10) - 5
        enddo
      enddo
      call sf(abcss, n_scatt, xyz, b_iso, n_refl, hkl, f_calc)
      if (n_refl .le. 100) then
        DO i=1,n_refl
          write(6, '(3(1x,i3),1x,f12.6,1x,f12.6)')
     &      hkl(1,i), hkl(2,i), hkl(3,i),
     &      f_calc(1,i), f_calc(2,i)
        enddo
      else
        max_a = 0
        max_b = 0
        DO i=1,n_refl
          a = f_calc(1,i)
          b = f_calc(2,i)
          if (max_a .lt. a) max_a = a
          if (max_b .lt. b) max_b = b
        enddo
        write(6, '(2(1x,f12.6))') max_a, max_b
      endif
      end
