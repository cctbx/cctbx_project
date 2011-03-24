#include <fem.hpp> // Fortran EMulation library of fable module

namespace cctbx {
namespace xray {
namespace targets {

using namespace fem::major_types;

using fem::common;

//C
//C  Differentiation of calc_k_dv in forward (tangent) mode: (multi-directional mode)
//C   variations   of useful results: k kd
//C   with respect to varying inputs: ic
//C
inline
void
calc_k_dv_dv(
  double& k,
  arr_ref<double> kd0,
  arr_ref<double> kd,
  arr_ref<double, 2> kdd,
  int const& nh,
  arr_cref<double> fo,
  arr_cref<double> ic,
  arr_cref<double, 2> icd0,
  arr_cref<double, 2> icd,
  int const& nbdirs,
  int const& nbdirs0)
{
  kd0(dimension(nbdirs));
  kd(dimension(nbdirs));
  kdd(dimension(nbdirs, nbdirs));
  fo(dimension(nh));
  ic(dimension(nh));
  icd0(dimension(nbdirs, nh));
  icd(dimension(nbdirs, nh));
  double k_num = 0;
  double k_den = 0;
  int nd = fem::int0;
  int nd0 = fem::int0;
  arr<double, 2> k_dendd(dimension(nbdirs, nbdirs), fem::fill0);
  arr<double> k_dend(dimension(nbdirs), fem::fill0);
  FEM_DO(nd, 1, nbdirs) {
    FEM_DO(nd0, 1, nbdirs0) {
      k_dendd(nd0, nd) = 0.e0;
    }
    k_dend(nd) = 0.e0;
  }
  arr<double, 2> k_numdd(dimension(nbdirs, nbdirs), fem::fill0);
  arr<double> k_numd(dimension(nbdirs), fem::fill0);
  FEM_DO(nd, 1, nbdirs) {
    FEM_DO(nd0, 1, nbdirs0) {
      k_numdd(nd0, nd) = 0.e0;
    }
    k_numd(nd) = 0.e0;
  }
  arr<double> k_dend0(dimension(nbdirs), fem::fill0);
  FEM_DO(nd0, 1, nbdirs0) {
    k_dend0(nd0) = 0.e0;
  }
  arr<double> k_numd0(dimension(nbdirs), fem::fill0);
  FEM_DO(nd0, 1, nbdirs0) {
    k_numd0(nd0) = 0.e0;
  }
  int ii1 = fem::int0;
  FEM_DO(nd0, 1, nbdirs0) {
    FEM_DO(ii1, 1, nbdirs) {
      k_numdd(nd0, ii1) = 0.e0;
    }
  }
  arr<double, 2> result1dd(dimension(nbdirs, nbdirs), fem::fill0);
  FEM_DO(nd0, 1, nbdirs0) {
    FEM_DO(ii1, 1, nbdirs) {
      result1dd(nd0, ii1) = 0.e0;
    }
  }
  int ih = fem::int0;
  arr<double> result1d(dimension(nbdirs), fem::fill0);
  double result10 = fem::double0;
  arr<double> result10d(dimension(nbdirs), fem::fill0);
  arr<double> result1d0(dimension(nbdirs), fem::fill0);
  double result1 = fem::double0;
  FEM_DO(ih, 1, nh) {
    FEM_DO(nd, 1, nbdirs) {
      if (ic(ih) == 0.0f) {
        FEM_DO(nd0, 1, nbdirs0) {
          result1dd(nd0, nd) = 0.e0;
        }
        result1d(nd) = 0.e0;
      }
      else {
        result10 = fem::sqrt(ic(ih));
        FEM_DO(nd0, 1, nbdirs0) {
          if (ic(ih) == 0.0f) {
            result10d(nd0) = 0.e0;
          }
          else {
            result10d(nd0) = icd0(nd0, ih) / (2.0f * fem::sqrt(ic(ih)));
          }
          result1dd(nd0, nd) = -(icd(nd, ih) * 2.0f * result10d(
            nd0) / fem::pow2((2.0f * result10)));
        }
        result1d(nd) = icd(nd, ih) / (2.0f * result10);
      }
      FEM_DO(nd0, 1, nbdirs0) {
        k_numdd(nd0, nd) += fo(ih) * result1dd(nd0, nd);
        k_dendd(nd0, nd) = 0.e0;
      }
      k_numd(nd) += fo(ih) * result1d(nd);
      k_dend(nd) += icd(nd, ih);
    }
    FEM_DO(nd0, 1, nbdirs0) {
      if (ic(ih) == 0.0f) {
        result1d0(nd0) = 0.e0;
      }
      else {
        result1d0(nd0) = icd0(nd0, ih) / (2.0f * fem::sqrt(ic(ih)));
      }
      k_numd0(nd0) += fo(ih) * result1d0(nd0);
      k_dend0(nd0) += icd0(nd0, ih);
    }
    result1 = fem::sqrt(ic(ih));
    k_num += fo(ih) * result1;
    k_den += ic(ih);
  }
  FEM_DO(nd0, 1, nbdirs0) {
    FEM_DO(ii1, 1, nbdirs) {
      kdd(nd0, ii1) = 0.e0;
    }
  }
  FEM_DO(nd, 1, nbdirs) {
    FEM_DO(nd0, 1, nbdirs0) {
      kdd(nd0, nd) = ((k_numdd(nd0, nd) * k_den + k_numd(nd) *
        k_dend0(nd0) - k_dend(nd) * k_numd0(nd0)) * fem::pow2(
        k_den) - (k_numd(nd) * k_den - k_num * k_dend(nd)) * 2 *
        k_den * k_dend0(nd0)) / fem::pow2((fem::pow2(k_den)));
    }
    kd(nd) = (k_numd(nd) * k_den - k_num * k_dend(nd)) / fem::pow2(k_den);
  }
  FEM_DO(nd0, 1, nbdirs0) {
    kd0(nd0) = (k_numd0(nd0) * k_den - k_num * k_dend0(nd0)) / fem::pow2(k_den);
  }
  k = k_num / k_den;
}

//C
//C  Differentiation of calc_w_dv in forward (tangent) mode: (multi-directional mode)
//C   variations   of useful results: w wd
//C   with respect to varying inputs: k kd ic
//C
inline
void
calc_w_dv_dv(
  arr_ref<double> w,
  arr_ref<double, 2> wd0,
  arr_ref<double, 2> wd,
  arr_ref<double, 3> wdd,
  int const& nh,
  arr_cref<double> io,
  arr_cref<double> so,
  arr_cref<double> ic,
  arr_cref<double, 2> icd0,
  arr_cref<double, 2> icd,
  double const& k,
  arr_cref<double> kd0,
  arr_cref<double> kd,
  arr_cref<double, 2> kdd,
  double const& wa,
  double const& wb,
  int const& nbdirs,
  int const& nbdirs0)
{
  w(dimension(nh));
  wd0(dimension(nbdirs, nh));
  wd(dimension(nbdirs, nh));
  wdd(dimension(nbdirs, nbdirs, nh));
  io(dimension(nh));
  so(dimension(nh));
  ic(dimension(nh));
  icd0(dimension(nbdirs, nh));
  icd(dimension(nbdirs, nh));
  kd0(dimension(nbdirs));
  kd(dimension(nbdirs));
  kdd(dimension(nbdirs, nbdirs));
  int nd0 = fem::int0;
  int ii10 = fem::int0;
  arr<double, 2> k_sqdd(dimension(nbdirs, nbdirs), fem::fill0);
  FEM_DO(nd0, 1, nbdirs0) {
    FEM_DO(ii10, 1, nbdirs) {
      k_sqdd(nd0, ii10) = 0.e0;
    }
  }
  int nd = fem::int0;
  arr<double> k_sqd(dimension(nbdirs), fem::fill0);
  FEM_DO(nd, 1, nbdirs) {
    FEM_DO(nd0, 1, nbdirs0) {
      k_sqdd(nd0, nd) = 2 * (kd0(nd0) * kd(nd) + k * kdd(nd0, nd));
    }
    k_sqd(nd) = 2 * k * kd(nd);
  }
  arr<double> k_sqd0(dimension(nbdirs), fem::fill0);
  FEM_DO(nd0, 1, nbdirs0) {
    k_sqd0(nd0) = 2 * k * kd0(nd0);
  }
  double k_sq = fem::pow2(k);
  int ii1 = fem::int0;
  FEM_DO(nd, 1, nbdirs) {
    FEM_DO(ii1, 1, nh) {
      FEM_DO(nd0, 1, nbdirs0) {
        wdd(nd0, nd, ii1) = 0.e0;
      }
      wd(nd, ii1) = 0.e0;
    }
  }
  FEM_DO(nd0, 1, nbdirs0) {
    FEM_DO(ii10, 1, nh) {
      wd0(nd0, ii10) = 0.e0;
    }
  }
  int ii2 = fem::int0;
  FEM_DO(nd0, 1, nbdirs0) {
    FEM_DO(ii10, 1, nh) {
      FEM_DO(ii2, 1, nbdirs) {
        wdd(nd0, ii2, ii10) = 0.e0;
      }
    }
  }
  arr<double, 2> ikdd(dimension(nbdirs, nbdirs), fem::fill0);
  FEM_DO(nd0, 1, nbdirs0) {
    FEM_DO(ii10, 1, nbdirs) {
      ikdd(nd0, ii10) = 0.e0;
    }
  }
  arr<double, 2> skdd(dimension(nbdirs, nbdirs), fem::fill0);
  FEM_DO(nd0, 1, nbdirs0) {
    FEM_DO(ii10, 1, nbdirs) {
      skdd(nd0, ii10) = 0.e0;
    }
  }
  arr<double, 2> wa_p_sqdd(dimension(nbdirs, nbdirs), fem::fill0);
  FEM_DO(nd0, 1, nbdirs0) {
    FEM_DO(ii10, 1, nbdirs) {
      wa_p_sqdd(nd0, ii10) = 0.e0;
    }
  }
  arr<double, 2> pdd(dimension(nbdirs, nbdirs), fem::fill0);
  FEM_DO(nd0, 1, nbdirs0) {
    FEM_DO(ii10, 1, nbdirs) {
      pdd(nd0, ii10) = 0.e0;
    }
  }
  arr<double, 2> sk_sqdd(dimension(nbdirs, nbdirs), fem::fill0);
  FEM_DO(nd0, 1, nbdirs0) {
    FEM_DO(ii10, 1, nbdirs) {
      sk_sqdd(nd0, ii10) = 0.e0;
    }
  }
  int ih = fem::int0;
  arr<double> ikd(dimension(nbdirs), fem::fill0);
  arr<double> skd(dimension(nbdirs), fem::fill0);
  arr<double> ikd0(dimension(nbdirs), fem::fill0);
  arr<double> skd0(dimension(nbdirs), fem::fill0);
  double ik = fem::double0;
  double sk = fem::double0;
  double p = fem::double0;
  arr<double> pd0(dimension(nbdirs), fem::fill0);
  arr<double> sk_sqd0(dimension(nbdirs), fem::fill0);
  arr<double> wa_p_sqd0(dimension(nbdirs), fem::fill0);
  double sk_sq = fem::double0;
  double wa_p_sq = fem::double0;
  arr<double> pd(dimension(nbdirs), fem::fill0);
  arr<double> sk_sqd(dimension(nbdirs), fem::fill0);
  arr<double> wa_p_sqd(dimension(nbdirs), fem::fill0);
  FEM_DO(ih, 1, nh) {
    FEM_DO(nd, 1, nbdirs) {
      FEM_DO(nd0, 1, nbdirs0) {
        ikdd(nd0, nd) = -((io(ih) * k_sqdd(nd0, nd) * fem::pow2(
          k_sq) - io(ih) * k_sqd(nd) * 2 * k_sq * k_sqd0(nd0)) /
          fem::pow2((fem::pow2(k_sq))));
        skdd(nd0, nd) = -((so(ih) * k_sqdd(nd0, nd) * fem::pow2(
          k_sq) - so(ih) * k_sqd(nd) * 2 * k_sq * k_sqd0(nd0)) /
          fem::pow2((fem::pow2(k_sq))));
      }
      ikd(nd) = -(io(ih) * k_sqd(nd) / fem::pow2(k_sq));
      skd(nd) = -(so(ih) * k_sqd(nd) / fem::pow2(k_sq));
    }
    FEM_DO(nd0, 1, nbdirs0) {
      ikd0(nd0) = -(io(ih) * k_sqd0(nd0) / fem::pow2(k_sq));
      skd0(nd0) = -(so(ih) * k_sqd0(nd0) / fem::pow2(k_sq));
    }
    ik = io(ih) / k_sq;
    sk = so(ih) / k_sq;
    if (ik < 0) {
      ik = 0;
      FEM_DO(nd, 1, nbdirs) {
        FEM_DO(nd0, 1, nbdirs0) {
          ikdd(nd0, nd) = 0.e0;
        }
        ikd(nd) = 0.e0;
      }
      FEM_DO(nd0, 1, nbdirs0) {
        ikd0(nd0) = 0.e0;
      }
    }
    p = (ik + 2 * ic(ih)) / 3;
    FEM_DO(nd0, 1, nbdirs0) {
      pd0(nd0) = (ikd0(nd0) + 2 * icd0(nd0, ih)) / 3;
      sk_sqd0(nd0) = 2 * sk * skd0(nd0);
      wa_p_sqd0(nd0) = 2 * fem::pow2(wa) * p * pd0(nd0);
    }
    sk_sq = fem::pow2(sk);
    wa_p_sq = fem::pow2((wa * p));
    FEM_DO(nd, 1, nbdirs) {
      pd(nd) = (ikd(nd) + 2 * icd(nd, ih)) / 3;
      sk_sqd(nd) = 2 * sk * skd(nd);
      wa_p_sqd(nd) = 2 * fem::pow2(wa) * p * pd(nd);
      FEM_DO(nd0, 1, nbdirs0) {
        pdd(nd0, nd) = ikdd(nd0, nd) / 3;
        sk_sqdd(nd0, nd) = 2 * (skd0(nd0) * skd(nd) + sk * skdd(nd0, nd));
        wa_p_sqdd(nd0, nd) = 2 * fem::pow2(wa) * (pd0(nd0) * pd(nd) +
          p * pdd(nd0, nd));
        wdd(nd0, nd, ih) = ((sk_sqd(nd) + wa_p_sqd(nd) + wb * pd(
          nd)) * 2 * (sk_sq + wa_p_sq + wb * p) * (sk_sqd0(nd0) +
          wa_p_sqd0(nd0) + wb * pd0(nd0)) - (sk_sqdd(nd0, nd) + wa_p_sqdd(nd0,
          nd) + wb * pdd(nd0, nd)) * fem::pow2((sk_sq + wa_p_sq +
          wb * p))) / fem::pow2((fem::pow2((sk_sq + wa_p_sq + wb *
          p))));
      }
      wd(nd, ih) = (-(sk_sqd(nd) + wa_p_sqd(nd) + wb * pd(nd))) /
        fem::pow2((sk_sq + wa_p_sq + wb * p));
    }
    FEM_DO(nd0, 1, nbdirs0) {
      wd0(nd0, ih) = (-(sk_sqd0(nd0) + wa_p_sqd0(nd0) + wb * pd0(
        nd0))) / fem::pow2((sk_sq + wa_p_sq + wb * p));
    }
    w(ih) = 1 / (sk_sq + wa_p_sq + wb * p);
  }
}

//C
//C  Differentiation of calc_t_dv in forward (tangent) mode: (multi-directional mode)
//C   variations   of useful results: t td
//C   with respect to varying inputs: k kd w ic wd
//C
inline
void
calc_t_dv_dv(
  double& t,
  arr_ref<double> td0,
  arr_ref<double> td,
  arr_ref<double, 2> tdd,
  int const& nh,
  arr_cref<double> io,
  arr_cref<double> ic,
  arr_cref<double, 2> icd0,
  arr_cref<double, 2> icd,
  double const& k,
  arr_cref<double> kd0,
  arr_cref<double> kd,
  arr_cref<double, 2> kdd,
  arr_cref<double> w,
  arr_cref<double, 2> wd0,
  arr_cref<double, 2> wd,
  arr_cref<double, 3> wdd,
  int const& nbdirs,
  int const& nbdirs0)
{
  td0(dimension(nbdirs));
  td(dimension(nbdirs));
  tdd(dimension(nbdirs, nbdirs));
  io(dimension(nh));
  ic(dimension(nh));
  icd0(dimension(nbdirs, nh));
  icd(dimension(nbdirs, nh));
  kd0(dimension(nbdirs));
  kd(dimension(nbdirs));
  kdd(dimension(nbdirs, nbdirs));
  w(dimension(nh));
  wd0(dimension(nbdirs, nh));
  wd(dimension(nbdirs, nh));
  wdd(dimension(nbdirs, nbdirs, nh));
  int nd0 = fem::int0;
  int ii1 = fem::int0;
  arr<double, 2> k_sqdd(dimension(nbdirs, nbdirs), fem::fill0);
  FEM_DO(nd0, 1, nbdirs0) {
    FEM_DO(ii1, 1, nbdirs) {
      k_sqdd(nd0, ii1) = 0.e0;
    }
  }
  int nd = fem::int0;
  arr<double> k_sqd(dimension(nbdirs), fem::fill0);
  FEM_DO(nd, 1, nbdirs) {
    FEM_DO(nd0, 1, nbdirs0) {
      k_sqdd(nd0, nd) = 2 * (kd0(nd0) * kd(nd) + k * kdd(nd0, nd));
    }
    k_sqd(nd) = 2 * k * kd(nd);
  }
  arr<double> k_sqd0(dimension(nbdirs), fem::fill0);
  FEM_DO(nd0, 1, nbdirs0) {
    k_sqd0(nd0) = 2 * k * kd0(nd0);
  }
  double k_sq = fem::pow2(k);
  double t_num = 0;
  double t_den = 0;
  arr<double, 2> t_dendd(dimension(nbdirs, nbdirs), fem::fill0);
  arr<double> t_dend(dimension(nbdirs), fem::fill0);
  FEM_DO(nd, 1, nbdirs) {
    FEM_DO(nd0, 1, nbdirs0) {
      t_dendd(nd0, nd) = 0.e0;
    }
    t_dend(nd) = 0.e0;
  }
  arr<double, 2> t_numdd(dimension(nbdirs, nbdirs), fem::fill0);
  arr<double> t_numd(dimension(nbdirs), fem::fill0);
  FEM_DO(nd, 1, nbdirs) {
    FEM_DO(nd0, 1, nbdirs0) {
      t_numdd(nd0, nd) = 0.e0;
    }
    t_numd(nd) = 0.e0;
  }
  arr<double> t_dend0(dimension(nbdirs), fem::fill0);
  FEM_DO(nd0, 1, nbdirs0) {
    t_dend0(nd0) = 0.e0;
  }
  arr<double> t_numd0(dimension(nbdirs), fem::fill0);
  FEM_DO(nd0, 1, nbdirs0) {
    t_numd0(nd0) = 0.e0;
  }
  FEM_DO(nd0, 1, nbdirs0) {
    FEM_DO(ii1, 1, nbdirs) {
      t_numdd(nd0, ii1) = 0.e0;
    }
  }
  FEM_DO(nd0, 1, nbdirs0) {
    FEM_DO(ii1, 1, nbdirs) {
      t_dendd(nd0, ii1) = 0.e0;
    }
  }
  int ih = fem::int0;
  FEM_DO(ih, 1, nh) {
    FEM_DO(nd, 1, nbdirs) {
      FEM_DO(nd0, 1, nbdirs0) {
        t_numdd(nd0, nd) += wdd(nd0, nd, ih) * fem::pow2((io(ih) -
          k_sq * ic(ih))) + wd(nd, ih) * 2 * (io(ih) - k_sq * ic(
          ih)) * (-(k_sqd0(nd0) * ic(ih)) - k_sq * icd0(nd0, ih)) +
          2 * ((wd0(nd0, ih) * (io(ih) - k_sq * ic(ih)) + w(ih) * (-(
          k_sqd0(nd0) * ic(ih)) - k_sq * icd0(nd0, ih))) * (-(k_sqd(
          nd) * ic(ih)) - k_sq * icd(nd, ih))) + 2 * (w(ih) * (io(
          ih) - k_sq * ic(ih)) * (-(k_sqdd(nd0, nd) * ic(ih)) - k_sqd(
          nd) * icd0(nd0, ih) - icd(nd, ih) * k_sqd0(nd0)));
        t_dendd(nd0, nd) += fem::pow2(io(ih)) * wdd(nd0, nd, ih);
      }
      t_numd(nd) += wd(nd, ih) * fem::pow2((io(ih) - k_sq * ic(
        ih))) + w(ih) * 2 * (io(ih) - k_sq * ic(ih)) * (-(k_sqd(nd) *
        ic(ih)) - k_sq * icd(nd, ih));
      t_dend(nd) += fem::pow2(io(ih)) * wd(nd, ih);
    }
    FEM_DO(nd0, 1, nbdirs0) {
      t_numd0(nd0) += wd0(nd0, ih) * fem::pow2((io(ih) - k_sq * ic(
        ih))) + w(ih) * 2 * (io(ih) - k_sq * ic(ih)) * (-(k_sqd0(
        nd0) * ic(ih)) - k_sq * icd0(nd0, ih));
      t_dend0(nd0) += fem::pow2(io(ih)) * wd0(nd0, ih);
    }
    t_num += w(ih) * fem::pow2((io(ih) - k_sq * ic(ih)));
    t_den += w(ih) * fem::pow2(io(ih));
  }
  FEM_DO(nd0, 1, nbdirs0) {
    FEM_DO(ii1, 1, nbdirs) {
      tdd(nd0, ii1) = 0.e0;
    }
  }
  FEM_DO(nd, 1, nbdirs) {
    FEM_DO(nd0, 1, nbdirs0) {
      tdd(nd0, nd) = ((t_numdd(nd0, nd) * t_den + t_numd(nd) *
        t_dend0(nd0) - t_numd0(nd0) * t_dend(nd) - t_num * t_dendd(nd0,
        nd)) * fem::pow2(t_den) - (t_numd(nd) * t_den - t_num *
        t_dend(nd)) * 2 * t_den * t_dend0(nd0)) / fem::pow2((
        fem::pow2(t_den)));
    }
    td(nd) = (t_numd(nd) * t_den - t_num * t_dend(nd)) / fem::pow2(t_den);
  }
  FEM_DO(nd0, 1, nbdirs0) {
    td0(nd0) = (t_numd0(nd0) * t_den - t_num * t_dend0(nd0)) / fem::pow2(t_den);
  }
  t = t_num / t_den;
}

//C        Generated by TAPENADE     (INRIA, Tropics team)
//C  Tapenade 3.5 (r3782) - 22 Mar 2011 14:14
//C
//C  Differentiation of kwt_dv in forward (tangent) mode: (multi-directional mode)
//C   variations   of useful results: t td
//C   with respect to varying inputs: ic
//C   RW status of diff variables: t:out ic:in td:out
inline
void
kwt_dv_dv(
  double& t,
  arr_ref<double> td0,
  arr_ref<double> td,
  arr_ref<double, 2> tdd,
  int const& nh,
  arr_cref<double> fo,
  arr_cref<double> io,
  arr_cref<double> so,
  arr_cref<double> ic,
  arr_cref<double, 2> icd0,
  arr_cref<double, 2> icd,
  double const& wa,
  double const& wb,
  int const& nbdirs,
  int const& nbdirs0)
{
  td0(dimension(nbdirs));
  td(dimension(nbdirs));
  tdd(dimension(nbdirs, nbdirs));
  fo(dimension(nh));
  io(dimension(nh));
  so(dimension(nh));
  ic(dimension(nh));
  icd0(dimension(nbdirs, nh));
  icd(dimension(nbdirs, nh));
  double k = fem::double0;
  arr<double> kd0(dimension(nbdirs), fem::fill0);
  arr<double> kd(dimension(nbdirs), fem::fill0);
  arr<double, 2> kdd(dimension(nbdirs, nbdirs), fem::fill0);
  calc_k_dv_dv(k, kd0, kd, kdd, nh, fo, ic, icd0, icd, nbdirs, nbdirs0);
  arr<double> w(dimension(nh), fem::fill0);
  arr<double, 2> wd0(dimension(nbdirs, nh), fem::fill0);
  arr<double, 2> wd(dimension(nbdirs, nh), fem::fill0);
  arr<double, 3> wdd(dimension(nbdirs, nbdirs, nh), fem::fill0);
  calc_w_dv_dv(w, wd0, wd, wdd, nh, io, so, ic, icd0, icd, k, kd0,
    kd, kdd, wa, wb, nbdirs, nbdirs0);
  calc_t_dv_dv(t, td0, td, tdd, nh, io, ic, icd0, icd, k, kd0, kd,
    kdd, w, wd0, wd, wdd, nbdirs, nbdirs0);
}

}}} // namespace cctbx::xray::targets
