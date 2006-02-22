#ifndef SCITBX_SCITBX_MATRIX_ROW_ECHELON_H
#define SCITBX_SCITBX_MATRIX_ROW_ECHELON_H

#include <scitbx/mat_ref.h>
#include <scitbx/array_family/small.h>
#include <scitbx/array_family/tiny.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/array_family/misc_functions.h>
#include <boost/rational.hpp>

namespace scitbx { namespace matrix { namespace row_echelon {

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
      SCITBX_ASSERT(t.begin() != 0 && t.n_rows() >= mr);
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
  back_substitution_int(
    scitbx::mat_const_ref<IntType> const& re_mx,
    const IntType* v = 0,
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
          scitbx::matrix::multiply(&re_mx(ir,icp), &sol[icp], 1, nv, 1,
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

  template <typename IntType, typename FloatType>
  bool
  back_substitution_float(
    scitbx::mat_const_ref<IntType> const& re_mx,
    const FloatType* v,
    FloatType* sol)
  {
    SCITBX_ASSERT(sol != 0);
    using std::size_t;
    size_t nr = re_mx.n_rows();
    size_t nc = re_mx.n_columns();
    for (size_t ir = nr; ir > 0;)
    {
      ir--;
      size_t ic;
      for(ic=0;ic<nc;ic++) {
        if (re_mx(ir,ic)) goto set_sol_ic;
      }
      if (v && v[ir] != 0) return false;
      continue;

      set_sol_ic:
                size_t icp = ic + 1;
      size_t nv = nc - icp;
      if (nv) {
        scitbx::matrix::multiply(&re_mx(ir,icp), &sol[icp], 1, nv, 1,
                                 &sol[ic]);
        sol[ic] = -sol[ic];
      }
      else {
        sol[ic] = 0;
      }
      if (v) sol[ic] += v[ir];
      sol[ic] /= static_cast<FloatType>(re_mx(ir,ic));
    }
    return true;
  }

  template <typename IntType, std::size_t MaxIndependentIndices = 6>
  struct independent
  {
    independent() {}

    independent(scitbx::mat_const_ref<IntType> const& re_mx)
    {
      using std::size_t;
      size_t nc = re_mx.n_columns();
      SCITBX_ASSERT(nc <= MaxIndependentIndices);
      bool flag_indep[MaxIndependentIndices];
      IntType* n_a = 0;
      SCITBX_ASSERT(back_substitution_int(re_mx, n_a, n_a, flag_indep) >= 1);
      for(size_t ic=0;ic<nc;ic++) {
        if (flag_indep[ic]) indices.push_back(ic);
      }
    }

    af::small<std::size_t, MaxIndependentIndices> indices;
  };

  template <typename NumType, unsigned MaxNRows, unsigned NCols>
  struct full_pivoting_small
  {
    af::tiny<NumType, NCols*NCols> echelon_form;
    af::small<unsigned, MaxNRows> row_perm;
    af::tiny<unsigned, NCols> col_perm;
    af::small<unsigned, NCols> pivot_cols;
    af::small<unsigned, NCols> free_cols;

    full_pivoting_small() {}

    full_pivoting_small(
      af::ref<NumType, af::c_grid<2> > const& m_work,
      NumType const& min_abs_pivot=0,
      unsigned max_rank=NCols)
    {
      SCITBX_ASSERT(m_work.accessor()[0] <= MaxNRows);
      SCITBX_ASSERT(m_work.accessor()[0] >= NCols);
      SCITBX_ASSERT(m_work.accessor()[1] == NCols);
      unsigned n_rows = m_work.accessor()[0];
      for(unsigned i=0;i<n_rows;i++) row_perm.push_back(i);
      for(unsigned i=0;i<NCols;i++) col_perm[i] = i;
      unsigned pr = 0;
      for(unsigned pc=0;pc<NCols;pc++) {
        // search for the next pivot value; here "m" is for "max"
        unsigned mr = pr;
        unsigned mc = pc;
        unsigned ir_nc = pr * NCols;
        NumType mv = m_work[ir_nc+pc];
        for(unsigned ir=pr;ir<n_rows;ir++,ir_nc+=NCols) {
          unsigned ir_ic = ir_nc + pc;
          for(unsigned ic=pc;ic<NCols;ic++) {
            NumType v = scitbx::fn::absolute(m_work[ir_ic++]);
            if (mv < v) {
              mv = v;
              mr = ir;
              mc = ic;
            }
          }
        }
        if (mv > min_abs_pivot && pr < max_rank) {
          if (mr != pr) swap_rows(m_work.begin(), pr, mr);
          if (mc != pc) swap_cols(m_work.begin(), pc, mc);
          // subtract multiple of pivot row from all rows below
          unsigned ir_nc = pr * NCols;
          unsigned pr_pc = ir_nc + pc;
          NumType v = m_work[pr_pc++];
          ir_nc += NCols;
          for(unsigned ir=pr+1;ir<n_rows;ir++,ir_nc+=NCols) {
            unsigned ir_ic = ir_nc + pc;
            NumType f = m_work[ir_ic] / v;
            m_work[ir_ic] = 0;
            ir_ic++;
            unsigned pr_ic = pr_pc;
            for(unsigned ic=pc+1;ic<NCols;ic++) {
              m_work[ir_ic++] -= f*m_work[pr_ic++];
            }
          }
          pr++;
          pivot_cols.push_back(pc);
        }
        else {
          free_cols.push_back(pc);
        }
        // copy result to local memory
        std::copy(
          m_work.begin(),
          m_work.begin()+(NCols*NCols),
          echelon_form.begin());
      }
    }

    af::tiny<NumType, NCols>
    back_substitution(af::small<NumType, NCols> const& free_values) const
    {
      SCITBX_ASSERT(free_values.size() == free_cols.size());
      af::tiny<NumType, NCols> perm_result;
      for(unsigned i=0;i<free_cols.size();i++) {
        perm_result[free_cols[i]] = free_values[i];
      }
      unsigned rank = pivot_cols.size();
      for(unsigned ip=0;ip<rank;ip++) {
        unsigned pr = rank-ip-1;
        unsigned pc = pivot_cols[pr];
        unsigned pr_pc = pr * NCols + pc;
        unsigned pr_ic = pr_pc + 1;
        NumType s = 0;
        for(unsigned ic=pc+1;ic<NCols;) {
          s -= perm_result[ic++] * echelon_form[pr_ic++];
        }
        perm_result[pc] = s / echelon_form[pr_pc];
      }
      af::tiny<NumType, NCols> result;
      for(unsigned i=0;i<NCols;i++) {
        result[col_perm[i]] = perm_result[i];
      }
      return result;
    }

    protected:

      void
      swap_rows(NumType* m_work, unsigned i, unsigned j)
      {
        unsigned ic = i*row_perm.size();
        unsigned jc = j*row_perm.size();
        for(unsigned c=0;c<NCols;c++) {
          std::swap(m_work[ic++], m_work[jc++]);
        }
        std::swap(row_perm[i], row_perm[j]);
      }

      void
      swap_cols(NumType* m_work, unsigned i, unsigned j)
      {
        unsigned ri = i;
        unsigned rj = j;
        for(unsigned r=0;r<row_perm.size();r++) {
          std::swap(m_work[ri], m_work[rj]);
          ri += NCols;
          rj += NCols;
        }
        std::swap(col_perm[i], col_perm[j]);
      }
  };

}}} // namespace scitbx::matrix::row_echelon

#endif // SCITBX_SCITBX_MATRIX_ROW_ECHELON_H
