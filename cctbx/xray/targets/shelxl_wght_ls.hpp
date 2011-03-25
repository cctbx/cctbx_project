#include <fem.hpp> // Fortran EMulation library of fable module

namespace cctbx {
namespace xray {
namespace targets {

using namespace fem::major_types;

using fem::common;

inline
void
calc_k(
  double& k,
  int const& nh,
  arr_cref<double> fo,
  arr_cref<double> ic)
{
  fo(dimension(nh));
  ic(dimension(nh));
  //C
  double k_num = 0;
  double k_den = 0;
  int ih = fem::int0;
  FEM_DO(ih, 1, nh) {
    k_num += fo(ih) * fem::sqrt(ic(ih));
    k_den += ic(ih);
  }
  k = k_num / k_den;
}

//C
//C-----
//C
inline
void
calc_w(
  arr_ref<double> w,
  int const& nh,
  arr_cref<double> io,
  arr_cref<double> so,
  arr_cref<double> ic,
  double const& k,
  double const& wa,
  double const& wb)
{
  w(dimension(nh));
  io(dimension(nh));
  so(dimension(nh));
  ic(dimension(nh));
  //C
  double k_sq = fem::pow2(k);
  int ih = fem::int0;
  double ik = fem::double0;
  double sk = fem::double0;
  double p = fem::double0;
  double sk_sq = fem::double0;
  double wa_p_sq = fem::double0;
  FEM_DO(ih, 1, nh) {
    ik = io(ih) / k_sq;
    sk = so(ih) / k_sq;
    if (ik < 0) {
      ik = 0;
    }
    p = (ik + 2 * ic(ih)) / 3;
    sk_sq = fem::pow2(sk);
    wa_p_sq = fem::pow2((wa * p));
    w(ih) = 1 / (sk_sq + wa_p_sq + wb * p);
  }
}

//C
//C-----
//C
inline
void
calc_t(
  double& t,
  int const& nh,
  arr_cref<double> io,
  arr_cref<double> ic,
  double const& k,
  arr_cref<double> w)
{
  io(dimension(nh));
  ic(dimension(nh));
  w(dimension(nh));
  //C
  double k_sq = fem::pow2(k);
  double t_num = 0;
  double t_den = 0;
  int ih = fem::int0;
  FEM_DO(ih, 1, nh) {
    t_num += w(ih) * fem::pow2((io(ih) - k_sq * ic(ih)));
    t_den += w(ih) * fem::pow2(io(ih));
  }
  t = t_num / t_den;
  //C
}

//C
//C-----
//C
inline
void
kwt(
  double& t,
  int const& nh,
  arr_cref<double> fo,
  arr_cref<double> io,
  arr_cref<double> so,
  arr_cref<double> ic,
  double const& wa,
  double const& wb)
{
  fo(dimension(nh));
  io(dimension(nh));
  so(dimension(nh));
  ic(dimension(nh));
  //C
  double k = fem::double0;
  calc_k(k, nh, fo, ic);
  arr<double> w(dimension(nh), fem::fill0);
  calc_w(w, nh, io, so, ic, k, wa, wb);
  calc_t(t, nh, io, ic, k, w);
  //C
}

//C
//C  Differentiation of calc_k in forward (tangent) mode: (multi-directional mode)
//C   variations   of useful results: k
//C   with respect to varying inputs: ic
inline
void
calc_k_dv(
  double& k,
  arr_ref<double> kd,
  int const& nh,
  arr_cref<double> fo,
  arr_cref<double> ic)
{
  int const nbdirs = nh;
  kd(dimension(nbdirs));
  fo(dimension(nh));
  ic(dimension(nh));
  //C
  double k_num = 0;
  double k_den = 0;
  int nd = fem::int0;
  arr<double> k_dend(dimension(nbdirs), fem::fill0);
  arr<double> k_numd(dimension(nbdirs), fem::fill0);
  int ih = fem::int0;
  arr<double> result1d(dimension(nbdirs), fem::fill0);
  double result1 = fem::double0;
  FEM_DO(ih, 1, nh) {
    nd = ih; {
      if (ic(ih) == 0.0f) {
        result1d(nd) = 0.e0;
      }
      else {
        result1d(nd) = 1 / (2.0f * fem::sqrt(ic(ih)));
      }
      k_numd(nd) += fo(ih) * result1d(nd);
      k_dend(nd) += 1;
    }
    result1 = fem::sqrt(ic(ih));
    k_num += fo(ih) * result1;
    k_den += ic(ih);
  }
  FEM_DO(nd, 1, nbdirs) {
    kd(nd) = (k_numd(nd) * k_den - k_num * k_dend(nd)) / fem::pow2(k_den);
  }
  k = k_num / k_den;
}

//C
//C  Differentiation of calc_w in forward (tangent) mode: (multi-directional mode)
//C   variations   of useful results: w
//C   with respect to varying inputs: k ic
//C
//C-----
//C
inline
void
calc_w_dv(
  arr_ref<double> w,
  arr_ref<double, 2> wd,
  int const& nh,
  arr_cref<double> io,
  arr_cref<double> so,
  arr_cref<double> ic,
  double const& k,
  arr_cref<double> kd,
  double const& wa,
  double const& wb)
{
  int const nbdirs = nh;
  w(dimension(nh));
  wd(dimension(nbdirs, nh));
  io(dimension(nh));
  so(dimension(nh));
  ic(dimension(nh));
  kd(dimension(nbdirs));
  //C
  int nd = fem::int0;
  arr<double> k_sqd(dimension(nbdirs), fem::fill0);
  FEM_DO(nd, 1, nbdirs) {
    k_sqd(nd) = 2 * k * kd(nd);
  }
  //C
  double k_sq = fem::pow2(k);
  int ii1 = fem::int0;
  FEM_DO(nd, 1, nbdirs) {
    ii1 = nd; {
      wd(nd, ii1) = 0.e0;
    }
  }
  int ih = fem::int0;
  arr<double> ikd(dimension(nbdirs), fem::fill0);
  arr<double> skd(dimension(nbdirs), fem::fill0);
  double ik = fem::double0;
  double sk = fem::double0;
  double p = fem::double0;
  double sk_sq = fem::double0;
  double wa_p_sq = fem::double0;
  arr<double> pd(dimension(nbdirs), fem::fill0);
  arr<double> sk_sqd(dimension(nbdirs), fem::fill0);
  arr<double> wa_p_sqd(dimension(nbdirs), fem::fill0);
  FEM_DO(ih, 1, nh) {
    FEM_DO(nd, 1, nbdirs) {
      ikd(nd) = -(io(ih) * k_sqd(nd) / fem::pow2(k_sq));
      skd(nd) = -(so(ih) * k_sqd(nd) / fem::pow2(k_sq));
    }
    ik = io(ih) / k_sq;
    sk = so(ih) / k_sq;
    if (ik < 0) {
      ik = 0;
      FEM_DO(nd, 1, nbdirs) {
        ikd(nd) = 0.e0;
      }
    }
    p = (ik + 2 * ic(ih)) / 3;
    sk_sq = fem::pow2(sk);
    wa_p_sq = fem::pow2((wa * p));
    FEM_DO(nd, 1, nbdirs) {
      pd(nd) = (ikd(nd) + (nd == ih ? 2 : 0)) / 3;
      sk_sqd(nd) = 2 * sk * skd(nd);
      wa_p_sqd(nd) = 2 * fem::pow2(wa) * p * pd(nd);
      wd(nd, ih) = (-(sk_sqd(nd) + wa_p_sqd(nd) + wb * pd(nd))) /
        fem::pow2((sk_sq + wa_p_sq + wb * p));
    }
    w(ih) = 1 / (sk_sq + wa_p_sq + wb * p);
  }
}

//C
//C  Differentiation of calc_k_b in forward (tangent) mode: (multi-directional mode)
//C   variations   of useful results: icb
//C   with respect to varying inputs: kb ic icb
//C
inline
void
calc_k_b_dv(
  double const& /* k */,
  double const& kb,
  arr_cref<double> kbd,
  int const& nh,
  arr_cref<double> fo,
  arr_cref<double> ic,
  arr_ref<double> icb,
  arr_ref<double, 2> icbd)
{
  int const nbdirs = nh;
  kbd(dimension(nbdirs));
  fo(dimension(nh));
  ic(dimension(nh));
  icb(dimension(nh));
  icbd(dimension(nbdirs, nh));
  double k_num = 0;
  double k_den = 0;
  int nd = fem::int0;
  arr<double> k_dend(dimension(nbdirs), fem::fill0);
  arr<double> k_numd(dimension(nbdirs), fem::fill0);
  int ih = fem::int0;
  arr<double> result1d(dimension(nbdirs), fem::fill0);
  double result1 = fem::double0;
  FEM_DO(ih, 1, nh) {
    nd = ih; {
      if (ic(ih) == 0.0f) {
        result1d(nd) = 0.e0;
      }
      else {
        result1d(nd) = 1 / (2.0f * fem::sqrt(ic(ih)));
      }
      k_numd(nd) += fo(ih) * result1d(nd);
      k_dend(nd) += 1;
    }
    result1 = fem::sqrt(ic(ih));
    k_num += fo(ih) * result1;
    k_den += ic(ih);
  }
  arr<double> k_numbd(dimension(nbdirs), fem::fill0);
  arr<double> k_denbd(dimension(nbdirs), fem::fill0);
  FEM_DO(nd, 1, nbdirs) {
    k_numbd(nd) = (kbd(nd) * k_den - kb * k_dend(nd)) / fem::pow2(k_den);
    k_denbd(nd) = -(((k_numd(nd) * kb + k_num * kbd(nd)) * fem::pow2(
      k_den) - k_num * kb * 2 * k_den * k_dend(nd)) / fem::pow2((
      fem::pow2(k_den))));
  }
  double k_numb = kb / k_den;
  double k_denb = -(k_num * kb / fem::pow2(k_den));
  FEM_DOSTEP(ih, nh, 1, -1) {
    if (ic(ih) == 0.0f) {
      nd = ih; {
        icbd(nd, ih) += k_denbd(nd);
      }
      icb(ih) += k_denb;
    }
    else {
      result1 = fem::sqrt(ic(ih));
      nd = ih; {
        if (ic(ih) == 0.0f) {
          result1d(nd) = 0.e0;
        }
        else {
          result1d(nd) = 1 / (2.0f * fem::sqrt(ic(ih)));
        }
        icbd(nd, ih) += (fo(ih) * k_numbd(nd) * 2.0f * result1 - fo(
          ih) * k_numb * 2.0f * result1d(nd)) / fem::pow2((2.0f *
          result1)) + k_denbd(nd);
      }
      icb(ih) += fo(ih) * k_numb / (2.0f * result1) + k_denb;
    }
  }
}

//C
//C  Differentiation of calc_w_b in forward (tangent) mode: (multi-directional mode)
//C   variations   of useful results: kb icb
//C   with respect to varying inputs: k kb wb0 ic icb
//C
inline
void
calc_w_b_dv(
  arr_cref<double> /* w */,
  arr_ref<double> wb0,
  arr_ref<double, 2> wb0d,
  int const& nh,
  arr_cref<double> io,
  arr_cref<double> so,
  arr_cref<double> ic,
  arr_ref<double> icb,
  arr_ref<double, 2> icbd,
  double const& k,
  arr_cref<double> kd,
  double& kb,
  arr_ref<double> kbd,
  double const& wa,
  double const& wb)
{
  int const nbdirs = nh;
  wb0(dimension(nh));
  wb0d(dimension(nbdirs, nh));
  io(dimension(nh));
  so(dimension(nh));
  ic(dimension(nh));
  icb(dimension(nh));
  icbd(dimension(nbdirs, nh));
  kd(dimension(nbdirs));
  kbd(dimension(nbdirs));
  int nd = fem::int0;
  arr<double> k_sqd(dimension(nbdirs), fem::fill0);
  FEM_DO(nd, 1, nbdirs) {
    k_sqd(nd) = 2 * k * kd(nd);
  }
  double k_sq = fem::pow2(k);
  double k_sqb = 0.e0;
  arr<double> k_sqbd(dimension(nbdirs), fem::fill0);
  int ih = fem::int0;
  arr<double> ikd(dimension(nbdirs), fem::fill0);
  double ik = fem::double0;
  double p = fem::double0;
  double sk = fem::double0;
  double sk_sq = fem::double0;
  double wa_p_sq = fem::double0;
  double temp = fem::double0;
  double tempb = fem::double0;
  double sk_sqb = fem::double0;
  double wa_p_sqb = fem::double0;
  arr<double> pd(dimension(nbdirs), fem::fill0);
  arr<double> skd(dimension(nbdirs), fem::fill0);
  arr<double> sk_sqd(dimension(nbdirs), fem::fill0);
  arr<double> wa_p_sqd(dimension(nbdirs), fem::fill0);
  arr<double> tempd(dimension(nbdirs), fem::fill0);
  arr<double> tempbd(dimension(nbdirs), fem::fill0);
  arr<double> sk_sqbd(dimension(nbdirs), fem::fill0);
  arr<double> wa_p_sqbd(dimension(nbdirs), fem::fill0);
  arr<double> pbd(dimension(nbdirs), fem::fill0);
  arr<double> skbd(dimension(nbdirs), fem::fill0);
  arr<double> ikbd(dimension(nbdirs), fem::fill0);
  double pb = fem::double0;
  double skb = fem::double0;
  double ikb = fem::double0;
  FEM_DOSTEP(ih, nh, 1, -1) {
    FEM_DO(nd, 1, nbdirs) {
      ikd(nd) = -(io(ih) * k_sqd(nd) / fem::pow2(k_sq));
    }
    ik = io(ih) / k_sq;
    if (ik < 0) {
      ik = 0;
      FEM_DO(nd, 1, nbdirs) {
        ikd(nd) = 0.e0;
      }
    }
    p = (ik + 2 * ic(ih)) / 3;
    sk = so(ih) / k_sq;
    sk_sq = fem::pow2(sk);
    wa_p_sq = fem::pow2((wa * p));
    temp = sk_sq + wa_p_sq + wb * p;
    tempb = -(wb0(ih) / fem::pow2(temp));
    sk_sqb = tempb;
    wa_p_sqb = tempb;
    FEM_DO(nd, 1, nbdirs) {
      pd(nd) = (ikd(nd) + (nd == ih ? 2 : 0)) / 3;
      skd(nd) = -(so(ih) * k_sqd(nd) / fem::pow2(k_sq));
      sk_sqd(nd) = 2 * sk * skd(nd);
      wa_p_sqd(nd) = 2 * fem::pow2(wa) * p * pd(nd);
      tempd(nd) = sk_sqd(nd) + wa_p_sqd(nd) + wb * pd(nd);
      tempbd(nd) = -((wb0d(nd, ih) * fem::pow2(temp) - wb0(ih) * 2 *
        temp * tempd(nd)) / fem::pow2((fem::pow2(temp))));
      sk_sqbd(nd) = tempbd(nd);
      wa_p_sqbd(nd) = tempbd(nd);
      pbd(nd) = 2 * fem::pow2(wa) * (pd(nd) * wa_p_sqb + p *
        wa_p_sqbd(nd)) + wb * tempbd(nd);
      wb0d(nd, ih) = 0.e0;
      skbd(nd) = 2 * (skd(nd) * sk_sqb + sk * sk_sqbd(nd));
      ikbd(nd) = pbd(nd) / 3;
      icbd(nd, ih) += 2 * pbd(nd) / 3;
    }
    pb = 2 * fem::pow2(wa) * p * wa_p_sqb + wb * tempb;
    wb0(ih) = 0.e0;
    skb = 2 * sk * sk_sqb;
    ikb = pb / 3;
    icb(ih) += 2 * pb / 3;
    if (ik < 0) {
      ikb = 0;
      FEM_DO(nd, 1, nbdirs) {
        ikbd(nd) = 0.e0;
      }
    }
    FEM_DO(nd, 1, nbdirs) {
      k_sqbd(nd) = k_sqbd(nd) - (io(ih) * ikbd(nd) * fem::pow2(
        k_sq) - io(ih) * ikb * 2 * k_sq * k_sqd(nd)) / fem::pow2((
        fem::pow2(k_sq))) - (so(ih) * skbd(nd) * fem::pow2(k_sq) - so(
        ih) * skb * 2 * k_sq * k_sqd(nd)) / fem::pow2((fem::pow2(
        k_sq)));
    }
    k_sqb = k_sqb - io(ih) * ikb / fem::pow2(k_sq) - so(ih) * skb /
      fem::pow2(k_sq);
  }
  FEM_DO(nd, 1, nbdirs) {
    kbd(nd) += 2 * (kd(nd) * k_sqb) + 2 * (k * k_sqbd(nd));
  }
  kb += 2 * k * k_sqb;
}

//C
//C  Differentiation of calc_t_b in forward (tangent) mode: (multi-directional mode)
//C   variations   of useful results: kb icb wb
//C   with respect to varying inputs: k w ic
//C
inline
void
calc_t_b_dv(
  double const& /* t */,
  double const& tb,
  int const& nh,
  arr_cref<double> io,
  arr_cref<double> ic,
  arr_ref<double> icb,
  arr_ref<double, 2> icbd,
  double const& k,
  arr_cref<double> kd,
  double& kb,
  arr_ref<double> kbd,
  arr_cref<double> w,
  arr_cref<double, 2> wd,
  arr_ref<double> wb,
  arr_ref<double, 2> wbd)
{
  int const nbdirs = nh;
  io(dimension(nh));
  ic(dimension(nh));
  icb(dimension(nh));
  icbd(dimension(nbdirs, nh));
  kd(dimension(nbdirs));
  kbd(dimension(nbdirs));
  w(dimension(nh));
  wd(dimension(nbdirs, nh));
  wb(dimension(nh));
  wbd(dimension(nbdirs, nh));
  int nd = fem::int0;
  arr<double> k_sqd(dimension(nbdirs), fem::fill0);
  FEM_DO(nd, 1, nbdirs) {
    k_sqd(nd) = 2 * k * kd(nd);
  }
  double k_sq = fem::pow2(k);
  double t_num = 0;
  double t_den = 0;
  arr<double> t_dend(dimension(nbdirs), fem::fill0);
  arr<double> t_numd(dimension(nbdirs), fem::fill0);
  int ih = fem::int0;
  FEM_DO(ih, 1, nh) {
    FEM_DO(nd, 1, nbdirs) {
      t_numd(nd) += wd(nd, ih) * fem::pow2((io(ih) - k_sq * ic(
        ih))) + w(ih) * 2 * (io(ih) - k_sq * ic(ih)) * (-(k_sqd(nd) *
        ic(ih)) - k_sq * (nd == ih ? 1 : 0));
      t_dend(nd) += fem::pow2(io(ih)) * wd(nd, ih);
    }
    t_num += w(ih) * fem::pow2((io(ih) - k_sq * ic(ih)));
    t_den += w(ih) * fem::pow2(io(ih));
  }
  arr<double> t_numbd(dimension(nbdirs), fem::fill0);
  arr<double> t_denbd(dimension(nbdirs), fem::fill0);
  FEM_DO(nd, 1, nbdirs) {
    t_numbd(nd) = -(tb * t_dend(nd) / fem::pow2(t_den));
    t_denbd(nd) = -((tb * t_numd(nd) * fem::pow2(t_den) - t_num *
      tb * 2 * t_den * t_dend(nd)) / fem::pow2((fem::pow2(t_den))));
  }
  double t_numb = tb / t_den;
  double t_denb = -(t_num * tb / fem::pow2(t_den));
  int ii1 = fem::int0;
  FEM_DO(ii1, 1, nh) {
    nd = ii1; {
      wbd(nd, ii1) = 0.e0;
    }
    wb(ii1) = 0.e0;
  }
  FEM_DO(ii1, 1, nh) {
    nd = ii1; {
      icbd(nd, ii1) = 0.e0;
    }
    icb(ii1) = 0.e0;
  }
  double k_sqb = 0.e0;
  int ii10 = fem::int0;
  FEM_DO(nd, 1, nbdirs) {
    ii10 = nd; {
      icbd(nd, ii10) = 0.e0;
    }
  }
  FEM_DO(nd, 1, nbdirs) {
    ii10 = nd; {
      wbd(nd, ii10) = 0.e0;
    }
  }
  arr<double> k_sqbd(dimension(nbdirs), fem::fill0);
  double temp = fem::double0;
  double tempb = fem::double0;
  arr<double> tempd(dimension(nbdirs), fem::fill0);
  arr<double> tempbd(dimension(nbdirs), fem::fill0);
  FEM_DOSTEP(ih, nh, 1, -1) {
    temp = io(ih) - k_sq * ic(ih);
    tempb = w(ih) * 2 * temp * t_numb;
    FEM_DO(nd, 1, nbdirs) {
      tempd(nd) = -(k_sqd(nd) * ic(ih)) - k_sq * (nd == ih ? 1 : 0);
      wbd(nd, ih) += 2 * temp * tempd(nd) * t_numb + fem::pow2(
        temp) * t_numbd(nd) + fem::pow2(io(ih)) * t_denbd(nd);
      tempbd(nd) = 2 * (wd(nd, ih) * temp * t_numb + w(ih) * (tempd(
        nd) * t_numb + temp * t_numbd(nd)));
      k_sqbd(nd) = k_sqbd(nd) - (nd == ih ? 1 : 0) * tempb - ic(ih)
        * tempbd(nd);
      icbd(nd, ih) = icbd(nd, ih) - k_sqd(nd) * tempb - k_sq * tempbd(nd);
    }
    wb(ih) += fem::pow2(temp) * t_numb + fem::pow2(io(ih)) * t_denb;
    k_sqb = k_sqb - ic(ih) * tempb;
    icb(ih) = icb(ih) - k_sq * tempb;
  }
  FEM_DO(nd, 1, nbdirs) {
    kbd(nd) = 2 * (kd(nd) * k_sqb + k * k_sqbd(nd));
  }
  kb = 2 * k * k_sqb;
}

//C        Generated by TAPENADE     (INRIA, Tropics team)
//C  Tapenade 3.5 (r3782) - 22 Mar 2011 14:14
//C
//C  Differentiation of kwt_b in forward (tangent) mode: (multi-directional mode)
//C   variations   of useful results: icb
//C   with respect to varying inputs: ic
//C   RW status of diff variables: ic:in icb:out
inline
void
kwt_b_dv(
  double const& t,
  double& tb,
  int const& nh,
  arr_cref<double> fo,
  arr_cref<double> io,
  arr_cref<double> so,
  arr_cref<double> ic,
  arr_ref<double> icb,
  arr_ref<double, 2> icbd,
  double const& wa,
  double const& wb)
{
  int const nbdirs = nh;
  fo(dimension(nh));
  io(dimension(nh));
  so(dimension(nh));
  ic(dimension(nh));
  icb(dimension(nh));
  icbd(dimension(nbdirs, nh));
  double k = fem::double0;
  arr<double> kd(dimension(nbdirs), fem::fill0);
  calc_k_dv(k, kd, nh, fo, ic);
  arr<double> w(dimension(nh), fem::fill0);
  arr<double, 2> wd(dimension(nbdirs, nh), fem::fill0);
  calc_w_dv(w, wd, nh, io, so, ic, k, kd, wa, wb);
  double kb = fem::double0;
  arr<double> kbd(dimension(nbdirs), fem::fill0);
  arr<double> wb0(dimension(nh), fem::fill0);
  arr<double, 2> wb0d(dimension(nbdirs, nh), fem::fill0);
  calc_t_b_dv(t, tb, nh, io, ic, icb, icbd, k, kd, kb, kbd, w,
    wd, wb0, wb0d);
  calc_w_b_dv(w, wb0, wb0d, nh, io, so, ic, icb, icbd, k, kd,
    kb, kbd, wa, wb);
  calc_k_b_dv(k, kb, kbd, nh, fo, ic, icb, icbd);
  tb = 0.e0;
}

}}} // namespace cctbx::xray::targets
