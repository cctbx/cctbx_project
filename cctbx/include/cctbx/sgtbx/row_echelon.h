/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Refactored parts of cctbx/math.h (rwgk)
     2001 May: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     2001 Apr: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_SGTBX_ROW_ECHELON_H
#define CCTBX_SGTBX_ROW_ECHELON_H

#include <cctbx/sgtbx/basic.h>
#include <scitbx/mat_ref.h>
#include <scitbx/array_family/small.h>
#include <scitbx/array_family/misc_functions.h>
#include <cctbx/import_scitbx_af.h>
#include <boost/rational.hpp>

namespace cctbx { namespace sgtbx { namespace row_echelon {

  namespace detail {

    template <typename AnyType>
    inline void
    swap(AnyType* a, AnyType* b, std::size_t n)
    {
      for(std::size_t i=0;i<n;i++) std::swap(a[i], b[i]);
    }

  } // namespace detail

  template <typename IntType>
  std::size_t
  form_t(scitbx::mat_ref<IntType>& m,
         scitbx::mat_ref<IntType> const& t)
  {
    // C++ version of RowEchelonFormT from the CrystGAP package
    //     (GAP Version 3.4.4).
    // B. Eick, F. Ga"hler and W. Nickel
    //     Computing Maximal Subgroups and Wyckoff Positions of Space Groups
    //     Acta Cryst. (1997). A53, 467 - 474
    using std::size_t;
    size_t mr = m.n_rows();
    size_t mc = m.n_columns();
    size_t tc = t.n_columns();
    if (tc) {
      CCTBX_ASSERT(t.begin() != 0 && t.n_rows() >= mr);
    }
    size_t i, j;
    for (i = j = 0; i < mr && j < mc;) {
      size_t k = i; while (k < mr && m(k,j) == 0) k++;
      if (k == mr)
        j++;
      else {
        if (i != k) {
                  detail::swap(&m(i,0), &m(k,0), mc);
          if (tc) detail::swap(&t(i,0), &t(k,0), tc);
        }
        for (k++; k < mr; k++) {
          IntType a = scitbx::fn::absolute(m(k, j));
          if (a != 0 && a < scitbx::fn::absolute(m(i,j))) {
                    detail::swap(&m(i,0), &m(k,0), mc);
            if (tc) detail::swap(&t(i,0), &t(k,0), tc);
          }
        }
        if (m(i,j) < 0) {
                  for(size_t ic=0;ic<mc;ic++) m(i,ic) *= -1;
          if (tc) for(size_t ic=0;ic<tc;ic++) t(i,ic) *= -1;
        }
        bool cleared = true;
        for (k = i+1; k < mr; k++) {
          IntType a = m(k,j) / m(i,j);
          if (a != 0) {
                    for(size_t ic=0;ic<mc;ic++) m(k,ic) -= a * m(i,ic);
            if (tc) for(size_t ic=0;ic<tc;ic++) t(k,ic) -= a * t(i,ic);
          }
          if (m(k,j) != 0) cleared = false;
        }
        if (cleared) { i++; j++; }
      }
    }
    m = scitbx::mat_ref<IntType>(m.begin(), i, mc);
    return i;
  }

  template <typename IntType>
  std::size_t
  form(scitbx::mat_ref<IntType>& m)
  {
    scitbx::mat_ref<IntType> t(0,0,0);
    return form_t(m, t);
  }

  template <typename IntType>
  IntType
  back_substitution(
    scitbx::mat_const_ref<IntType> const& re_mx,
    IntType* v = 0,
    IntType* sol = 0,
    bool* flag_indep = 0)
  {
    using std::size_t;
    size_t nr = re_mx.n_rows();
    size_t nc = re_mx.n_columns();
    if (flag_indep) {
      for(size_t ic=0;ic<nc;ic++) flag_indep[ic] = true;
    }
    IntType d = 1;
    for (size_t ir = nr; ir > 0;)
    {
      ir--;
      size_t ic;
      for(ic=0;ic<nc;ic++) {
        if (re_mx(ir,ic)) goto set_sol_ic;
      }
      if (v && v[ir] != 0) return 0;
      continue;

      set_sol_ic:
      if (flag_indep) flag_indep[ic] = false;
      if (sol) {
                  size_t icp = ic + 1;
        size_t nv = nc - icp;
        if (nv) {
          scitbx::matrix_multiply(&re_mx(ir,icp), &sol[icp], 1, nv, 1,
                                  &sol[ic]);
          sol[ic] *= -1;
        }
        else {
          sol[ic] = 0;
        }
        if (v) sol[ic] += d * v[ir];
                         IntType mrc = re_mx(ir,ic);
        IntType f = boost::gcd(sol[ic], mrc);
        if (mrc < 0) f *= -1;
        sol[ic] /= f;
        f = mrc / f;
        if (f != 1) {
          for(size_t jc=0;jc<nc;jc++) if (jc != ic) sol[jc] *= f;
          d *= f;
        }
      }
    }
    return d;
  }

  template <typename IntType, std::size_t MaxIndependentIndices = 6>
  struct independent
  {
    independent() {}

    independent(scitbx::mat_const_ref<IntType> const& re_mx)
    {
      using std::size_t;
      size_t nc = re_mx.n_columns();
      CCTBX_ASSERT(nc <= MaxIndependentIndices);
      bool flag_indep[MaxIndependentIndices];
      IntType* n_a = 0;
      CCTBX_ASSERT(back_substitution(re_mx, n_a, n_a, flag_indep) >= 1);
      for(size_t ic=0;ic<nc;ic++) {
        if (flag_indep[ic]) indices.push_back(ic);
      }
    }

    af::small<std::size_t, MaxIndependentIndices> indices;
  };

}}} // namespace cctbx::sgtbx::row_echelon

#endif // CCTBX_SGTBX_ROW_ECHELON_H
