#ifndef SCITBX_MATRIX_ROW_ECHELON_FULL_PIVOTING_IMPL_H
#define SCITBX_MATRIX_ROW_ECHELON_FULL_PIVOTING_IMPL_H

#include <algorithm>

namespace scitbx { namespace matrix { namespace row_echelon {
namespace full_pivoting_impl
{
  template <typename NumType>
  void
  swap_rows(
    NumType* a_work,
    unsigned n_cols,
    unsigned i,
    unsigned j)
  {
    unsigned ic = i*n_cols;
    unsigned jc = j*n_cols;
    for(unsigned c=0;c<n_cols;c++) {
      std::swap(a_work[ic++], a_work[jc++]);
    }
  }

  template <typename NumType>
  void
  swap_cols(
    NumType* a_work,
    unsigned n_rows,
    unsigned n_cols,
    unsigned i,
    unsigned j,
    unsigned* col_perm)
  {
    unsigned ri = i;
    unsigned rj = j;
    for(unsigned r=0;r<n_rows;r++) {
      std::swap(a_work[ri], a_work[rj]);
      ri += n_cols;
      rj += n_cols;
    }
    std::swap(col_perm[i], col_perm[j]);
  }

  //! For solving a * x = b, permitting n_rows != n_cols.
  template <typename NumType>
  unsigned
  reduction(
    unsigned n_rows,
    unsigned n_cols,
    NumType* a_work,
    NumType* b_work,
    NumType const& min_abs_pivot,
    unsigned max_rank,
    unsigned* col_perm) // n_cols
  {
    for(unsigned i=0;i<n_cols;i++) col_perm[i] = i;
    unsigned min_n_cols_n_rows = std::min(n_cols, n_rows);
    unsigned pr = 0;
    unsigned pc = 0;
    for(;pc<min_n_cols_n_rows;pc++) {
      // search for the next pivot value; here "m" is for "max"
      unsigned mr = pr;
      unsigned mc = pc;
      unsigned ir_nc = pr * n_cols;
      NumType mv = a_work[ir_nc+pc];
      for(unsigned ir=pr;ir<n_rows;ir++,ir_nc+=n_cols) {
        unsigned ir_ic = ir_nc + pc;
        for(unsigned ic=pc;ic<n_cols;ic++) {
          NumType v = a_work[ir_ic++];
          if (v < 0) v = -v;
          if (mv < v) {
            mv = v;
            mr = ir;
            mc = ic;
          }
        }
      }
      if (mv > min_abs_pivot && pr < max_rank) {
        if (mr != pr) {
          swap_rows(a_work, n_cols, pr, mr);
          if (b_work != 0) std::swap(b_work[pr], b_work[mr]);
        }
        if (mc != pc) {
          swap_cols(a_work, n_rows, n_cols, pc, mc, col_perm);
        }
        // subtract multiple of pivot row from all rows below
        unsigned ir_nc = pr * n_cols;
        unsigned pr_pc = ir_nc + pc;
        NumType v = a_work[pr_pc++];
        ir_nc += n_cols;
        for(unsigned ir=pr+1;ir<n_rows;ir++,ir_nc+=n_cols) {
          unsigned ir_ic = ir_nc + pc;
          NumType f = a_work[ir_ic] / v;
          a_work[ir_ic] = 0;
          ir_ic++;
          unsigned pr_ic = pr_pc;
          for(unsigned ic=pc+1;ic<n_cols;ic++) {
            a_work[ir_ic++] -= f*a_work[pr_ic++];
          }
          if (b_work != 0) b_work[ir] -= f*b_work[pr];
        }
        pr++;
      }
      else {
        break;
      }
    }
    return pr;
  }

  template <typename NumType>
  bool
  is_in_row_space(
    unsigned n_cols,
    const NumType* a_work,
    const unsigned* col_perm, // n_cols
    unsigned rank,
    NumType* x,
    NumType const& epsilon)
  {
    unsigned pr_pc = 0;
    for (unsigned i=0; i < rank; i++) {
      NumType xa = x[col_perm[i]] / a_work[pr_pc];
      if (xa != 0) {
        x[col_perm[i]] = 0;
        pr_pc++;
        for (unsigned k=i+1; k < n_cols; k++) {
          x[col_perm[k]] -= xa * a_work[pr_pc++];
        }
        pr_pc += i+1;
      }
      else {
        pr_pc += n_cols+1;
      }
    }
    for (unsigned i=0; i < n_cols; i++) {
      if (x[i] > epsilon) return false;
      if (x[i] < -epsilon) return false;
    }
    return true;
  }

  template <typename NumType>
  bool
  back_substitution(
    unsigned n_rows,
    unsigned n_cols,
    const NumType* a_work,
    const NumType* b_work,
    const unsigned* col_perm, // n_cols
    unsigned rank,
    const NumType* free_values,
    const NumType& epsilon,
    NumType* perm_result, // n_cols, scratch space
    NumType* result) // n_cols
  {
    if (b_work != 0) {
      for(unsigned i=rank;i<n_rows;i++) {
        if (b_work[i] > epsilon) return false;
        if (b_work[i] < -epsilon) return false;
      }
    }
    for(unsigned i=0;i<n_cols-rank;i++) {
      perm_result[rank+i] = free_values[i];
    }
    for(unsigned ip=0;ip<rank;ip++) {
      unsigned pr = rank-ip-1;
      unsigned pc = pr; // just for clarity
      unsigned pr_pc = pr * n_cols + pc;
      unsigned pr_ic = pr_pc + 1;
      NumType s = (b_work == 0 ? 0 : b_work[pr]);
      for(unsigned ic=pc+1;ic<n_cols;) {
        s -= perm_result[ic++] * a_work[pr_ic++];
      }
      perm_result[pc] = s / a_work[pr_pc];
    }
    for(unsigned i=0;i<n_cols;i++) {
      result[col_perm[i]] = perm_result[i];
    }
    return true;
  }

}}}} // namespace scitbx::matrix::row_echelon::full_pivoting_impl

#endif // include guard
