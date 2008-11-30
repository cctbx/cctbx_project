#ifndef SCITBX_MATRIX_ROW_ECHELON_FULL_PIVOTING_IMPL_H
#define SCITBX_MATRIX_ROW_ECHELON_FULL_PIVOTING_IMPL_H

#include <algorithm>

namespace scitbx { namespace matrix { namespace row_echelon {
namespace full_pivoting_impl
{
  template <typename NumType>
  void
  swap_rows(
    NumType* m_work,
    unsigned n_cols,
    unsigned i,
    unsigned j,
    unsigned* row_perm)
  {
    unsigned ic = i*n_cols;
    unsigned jc = j*n_cols;
    for(unsigned c=0;c<n_cols;c++) {
      std::swap(m_work[ic++], m_work[jc++]);
    }
    std::swap(row_perm[i], row_perm[j]);
  }

  template <typename NumType>
  void
  swap_cols(
    NumType* m_work,
    unsigned n_rows,
    unsigned n_cols,
    unsigned i,
    unsigned j,
    unsigned* col_perm)
  {
    unsigned ri = i;
    unsigned rj = j;
    for(unsigned r=0;r<n_rows;r++) {
      std::swap(m_work[ri], m_work[rj]);
      ri += n_cols;
      rj += n_cols;
    }
    std::swap(col_perm[i], col_perm[j]);
  }

  template <typename NumType>
  unsigned
  reduction(
    NumType* m_work,
    unsigned n_rows,
    unsigned n_cols,
    NumType const& min_abs_pivot,
    unsigned max_rank,
    unsigned* row_perm, // n_rows
    unsigned* col_perm, // n_cols
    unsigned* pivot_cols, // n_cols
    unsigned* free_cols) // n_cols
  {
    for(unsigned i=0;i<n_rows;i++) row_perm[i] = i;
    for(unsigned i=0;i<n_cols;i++) col_perm[i] = i;
    unsigned min_n_cols_n_rows = std::min(n_cols, n_rows);
    unsigned pr = 0;
    for(unsigned pc=0;pc<min_n_cols_n_rows;pc++) {
      // search for the next pivot value; here "m" is for "max"
      unsigned mr = pr;
      unsigned mc = pc;
      unsigned ir_nc = pr * n_cols;
      NumType mv = m_work[ir_nc+pc];
      for(unsigned ir=pr;ir<n_rows;ir++,ir_nc+=n_cols) {
        unsigned ir_ic = ir_nc + pc;
        for(unsigned ic=pc;ic<n_cols;ic++) {
          NumType v = m_work[ir_ic++];
          if (v < 0) v = -v;
          if (mv < v) {
            mv = v;
            mr = ir;
            mc = ic;
          }
        }
      }
      if (mv > min_abs_pivot && pr < max_rank) {
        if (mr != pr) swap_rows(m_work, n_cols, pr, mr, row_perm);
        if (mc != pc) swap_cols(m_work, n_rows, n_cols, pc, mc, col_perm);
        // subtract multiple of pivot row from all rows below
        unsigned ir_nc = pr * n_cols;
        unsigned pr_pc = ir_nc + pc;
        NumType v = m_work[pr_pc++];
        ir_nc += n_cols;
        for(unsigned ir=pr+1;ir<n_rows;ir++,ir_nc+=n_cols) {
          unsigned ir_ic = ir_nc + pc;
          NumType f = m_work[ir_ic] / v;
          m_work[ir_ic] = 0;
          ir_ic++;
          unsigned pr_ic = pr_pc;
          for(unsigned ic=pc+1;ic<n_cols;ic++) {
            m_work[ir_ic++] -= f*m_work[pr_ic++];
          }
        }
        pivot_cols[pr] = pc;
        pr++;
      }
      else {
        free_cols[pc-pr] = pc;
      }
    }
    for(unsigned pc=min_n_cols_n_rows;pc<n_cols;pc++) {
      free_cols[pc-pr] = pc;
    }
    return pr;
  }

  template <typename NumType>
  bool
  is_in_row_span(
    unsigned n_cols,
    const NumType* echelon_form,
    const unsigned* col_perm, // n_cols
    const unsigned* pivot_cols,
    unsigned pivot_cols_size,
    NumType* vector,
    NumType const& epsilon)
  {
    unsigned pr_pc = 0;
    for (unsigned i=0; i < pivot_cols_size; i++) {
      NumType a = vector[col_perm[i]] / echelon_form[pr_pc];
      if (a != 0) {
        for (int k=i; k < n_cols; k++) {
          vector[col_perm[k]] -= a * echelon_form[pr_pc++];
        }
      }
      pr_pc += i+1;
    }
    for (unsigned i=0; i < n_cols; i++) {
      if (vector[i] > epsilon) return false;
      if (vector[i] < -epsilon) return false;
    }
    return true;
  }

  template <typename NumType>
  void
  back_substitution(
    unsigned n_cols,
    const NumType* echelon_form,
    const unsigned* col_perm, // n_cols
    const unsigned* pivot_cols,
    unsigned pivot_cols_size,
    const unsigned* free_cols,
    unsigned free_cols_size,
    const NumType* free_values,
    NumType* perm_result, // n_cols, scratch space
    NumType* result) // n_cols
  {
    for(unsigned i=0;i<free_cols_size;i++) {
      perm_result[free_cols[i]] = free_values[i];
    }
    for(unsigned ip=0;ip<pivot_cols_size;ip++) {
      unsigned pr = pivot_cols_size-ip-1;
      unsigned pc = pivot_cols[pr];
      unsigned pr_pc = pr * n_cols + pc;
      unsigned pr_ic = pr_pc + 1;
      NumType s = 0;
      for(unsigned ic=pc+1;ic<n_cols;) {
        s -= perm_result[ic++] * echelon_form[pr_ic++];
      }
      perm_result[pc] = s / echelon_form[pr_pc];
    }
    for(unsigned i=0;i<n_cols;i++) {
      result[col_perm[i]] = perm_result[i];
    }
  }

}}}} // namespace scitbx::matrix::row_echelon::full_pivoting_impl

#endif // include guard
