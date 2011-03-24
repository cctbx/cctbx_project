      subroutine calc_k(k, nh, fo, ic)
      implicit none
c
      double precision k
      integer nh
      double precision fo(nh), ic(nh)
c
      integer ih
      double precision k_num, k_den
c
      k_num = 0
      k_den = 0
      do ih=1,nh
        k_num = k_num + fo(ih)*sqrt(ic(ih))
        k_den = k_den + ic(ih)
      enddo
      k = k_num/k_den
      end
c
c-----
c
      subroutine calc_w(w, nh, io, so, ic, k, wa, wb)
      implicit none
c
      double precision w(nh)
      integer nh
      double precision io(nh), so(nh), ic(nh), k, wa, wb
c
      integer ih
      double precision k_sq, ik, sk, p, sk_sq, wa_p_sq
c
      k_sq = k**2
      do ih=1,nh
        ik = io(ih)/k_sq
        sk = so(ih)/k_sq
        if (ik .lt. 0) ik = 0
        p = (ik+2*ic(ih))/3
        sk_sq = sk**2
        wa_p_sq = (wa*p)**2
        w(ih) = 1/(sk_sq+wa_p_sq+wb*p)
      enddo
      end
c
c-----
c
      subroutine calc_t(t, nh, io, ic, k, w)
      implicit none
c
      double precision t
      integer nh
      double precision io(nh), ic(nh), k, w(nh)
c
      integer ih
      double precision t_num, t_den, k_sq
c
      k_sq = k**2
      t_num = 0
      t_den = 0
      do ih=1,nh
        t_num = t_num + w(ih)*(io(ih)-k_sq*ic(ih))**2
        t_den = t_den + w(ih)*io(ih)**2
      enddo
      t = t_num / t_den
c
      end
c
c-----
c
      subroutine kwt(t, nh, fo, io, so, ic, wa, wb)
      implicit none
c
      double precision t
      integer nh
      double precision fo(nh), io(nh), so(nh), ic(nh), wa, wb
c
      double precision k
      double precision w(nh)
c
      call calc_k(k, nh, fo, ic)
      call calc_w(w, nh, io, so, ic, k, wa, wb)
      call calc_t(t, nh, io, ic, k, w)
c
      end
