#ifndef SCITBX_SPARSE_LU_DECOMPOSE_H
#define SCITBX_SPARSE_LU_DECOMPOSE_H

#include <cmath>
#include <iterator>
#include <algorithm>

#include <scitbx/sparse/dfs.h>

namespace scitbx { namespace sparse {

//! LU decomposition with partial pivoting.
/*!
 The algorithm is that of Gilbert and Peierls [1], [2, section 3.2.1].
 Its strength resides in running in a time proportional to the number of
 arithmetic operations involving non-zero elements, i.e. that the time spent
 working out the sparsity structure of L and U is kept proportional to the
 time necessary for those arithmetic operations (a remarkable property).

 [1] J.R. Gilbert and T. Peierls. Sparse partial pivoting in time
 proportional to arithmetic operations.
 SIAM Journal on Scientific and Statistical Computing, 9:862–874, 1988.

 [2] I.S. Duff and J.K. Reid. The design of ma48, a code for the direct
 solution of sparse unsymmetric linear systems of equations.
 Technical Report RAL-TR-95-039,
 Computer and Information Systems Department,
 Rutherford Appleton Laboratory, Oxon, UK, 1995.
*/
template<class Matrix>
class gilbert_peierls_lu_factorization
{
  public:
    typedef Matrix matrix_type;
    typedef typename Matrix::value_type value_type;
    typedef typename Matrix::column_type column_type;
    typedef typename Matrix::index_type index_type;
    typedef typename Matrix::row_iterator row_iterator;

    /// Construct the LU factorisation of m
    /** The columns of m may end up being permuted */
    gilbert_peierls_lu_factorization(const Matrix& m);

    /// The matrix this has been constructed with
    Matrix const& factored() {
      return a;
    }

    Matrix const& l() {
      return L;
    }

    Matrix const& u() {
      return U;
    }

    af::shared<index_type> rows_permutation() {
      return p;
    }

  private:
    /* We will use the notations of Golub and Van Loan (3rd edition)
    instead of those of [1],
    except that in accordance with C++ usage we use the notation 0:r to
    mean 1:r-1 in the Golub and Van Loan.
    Of special interest is algo (3.4.2) section 3.4.5.
    Modifications: we have m rows which may be different from the number n of
    columns; instead of swapping the rows physically, we will instead
    use indirection through the permutation.
    */

    /* We seek PA = LU */
    const Matrix& a;
    Matrix L, U;
    af::shared<index_type> p; /* Record of P:
                                    P(r) = s means that row r of A is
                                    row s of PA */
    af::shared<index_type> p_inv; /* Record of P^-1:
                                        P_inv(s) = r means that row r of A is
                                        row s of PA */

    typedef typename std::vector<index_type>::reverse_iterator
            row_idx_reverse_iter;
    typedef typename std::vector<index_type>::iterator
            row_idx_iter;

    /* Solve L(0:j, 0:j) z = A(0:j, j) for z,
      This is a forward substitution (c.f. Golub and Van Loan Algorithm 3.1.3)
      Then compute v(j:m) = A(j:m, j) - L(j:m, 0:j) z
    */
    // Symbolic phase, searching the sparsity
    void compute_z_and_v_sparsity(index_type j);
    std::vector<index_type> z_nz; // non-zero elements in z = U(0:j, j)
    std::vector<index_type> v_nz; // non-zero elements in v(j:m)

    // Numerical phase
    void compute_z_and_v(index_type j);

    // The DFS used to compute sparsity
    depth_first_search<Matrix> dfs;

    // The visitor for the computation of the sparsity of w
    struct w_sparsity;

    /* This starts as  w = 0 at each iteration over the columns
    (this is a dense vector but we will never touch those elements we can
    predict to be zero)
    then w contains the nonzero of A(:,j)
    then w = [ z ]
             [ v ]
    and finally
      w = [ U(0:j+1, j) ]
          [ L(j+1:m, j) ]
      (with the obvious clipping if j > m)
    */
    af::shared<value_type> w;

    // Copy the non-zero elements of A(:,j) into w = 0
    void initialize_w(index_type j);

    // Pivoting
    index_type find_pivot(index_type j);
    void swap_rows(index_type i, index_type j);

    // Write w into U(:,j) and L(:,j)
    void copy_z_into_U(index_type j);
    void copy_v_into_L(index_type j);
};


template<class Matrix>
gilbert_peierls_lu_factorization<Matrix>::
gilbert_peierls_lu_factorization(const Matrix& m)
  : a(m),
    L(a.n_rows(), std::min(a.n_rows(), a.n_cols())),
    U(std::min(a.n_rows(), a.n_cols()), a.n_cols()),
    p(m.n_rows()), p_inv(m.n_rows()),
    dfs(m.n_rows(), m.n_cols()),
    w(m.n_rows(), 0)
{
  // Initialise P to the identical permuntation
  for (index_type i=0; i < a.n_rows(); i++) p_inv[i] = p[i] = i;

  for (index_type j=0; j < U.n_cols(); j++) {
    // Computing column j of L and U
    a.col(j).fill_dense_vector_with_permutation(w, p);
    compute_z_and_v_sparsity(j);
    compute_z_and_v(j);
    if (j < L.n_rows() - 1) {
      index_type mu = find_pivot(j);
      swap_rows(mu, j);
    }
    copy_z_into_U(j);
    if (j < U.n_rows()) U(j,j) = w[j];
    if (j < L.n_rows()) copy_v_into_L(j);
  }

  // The rows of L have been kept unpermuted. Time to address that.
  L.permute_rows(p);
}


// Proxy to key data with a visitor interface to be plugged into a DFS
template<class Matrix>
struct gilbert_peierls_lu_factorization<Matrix>::w_sparsity
{
  index_type j;
  af::shared<index_type> &p;
  std::vector<index_type> &z_nz, &v_nz;

  w_sparsity(index_type j_,
             af::shared<index_type> &p_,
             std::vector<index_type> &z_nz_,
             std::vector<index_type>  &v_nz_)
    : j(j_), p(p_), z_nz(z_nz_), v_nz(v_nz_)
  {}

  index_type permute_rhs(index_type r) { return p[r]; }

  index_type permute(index_type r) { return permute_rhs(r); }

  // Empty the vectors which will hold the results on DFS start
  void dfs_started() {
    z_nz.clear();
    v_nz.clear();
  }

  /* A nonzero v(k) results from a nonzero A(k,j) */
  void dfs_started_from_vertex(index_type k) {
    if (k >= j) v_nz.push_back(k);
  }

  /* Only L(0:j, 0:j), so don't start a DFS from a z(k) with k >= j */
  bool dfs_shall_cut_tree_rooted_at(index_type k) {
    return k >= j;
  }

  /* A nonzero v(k) results from a nonzero z(l) through a nonzero L(k,l)
  for k >= j (the edge is l --> k, i.e. from column to row)
  */
  void dfs_found_tree_edge(index_type l, index_type k) {
    if (k >= j) v_nz.push_back(k);
  }

  /* A nonzero z(k) results from a nonzero z(l) through a nonzero L(k,l)
  for k < j and k != l.
  Hence don't continue the DFS further is k >= j or k == l.
  */
  bool dfs_shall_cut_tree_edge(index_type l, index_type k) {
    return k >= j || k == l;
  }

  /* All descendant of l have been found: thus z(k) does not depend on
  the z(p) for the indices p already in z_nz (but those z(p) depends on z(k))
  Hence a reverse iteration will give the indices in a topological order
  which makes the forward substitution to work.
  */
  void dfs_finished_vertex(index_type k) {
    z_nz.push_back(k);
  }
};

template<class Matrix>
inline
void gilbert_peierls_lu_factorization<Matrix>::
compute_z_and_v_sparsity(index_type j)
{

  /* Depth-first search
  Computing the sparsity of v and z at the same time is a trick from [2].
  */
  w_sparsity vis(j, p, z_nz, v_nz);
  dfs(L, a.col(j), vis);
}

template<class Matrix>
inline
void gilbert_peierls_lu_factorization<Matrix>::
compute_z_and_v(index_type j)
{
  // Forward substitution to compute z = w(0:j)
  for (row_idx_reverse_iter l_ = z_nz.rbegin(); l_ != z_nz.rend(); l_++) {
    index_type l = *l_;
    for (row_iterator k_ = L.col(l).begin(); k_ != L.col(l).end(); k_++) {
      index_type k = p[ k_.index() ]; // indirection
      if (k > j-1) continue; // working with L(0:j,0:j)
      if (k <= l) continue; // the j+1:n in Golub & Van Loan algo 3.1.3
      value_type L_kl = *k_;
      w[k] -= w[l] * L_kl;
    }
  }

  // Compute v(j:m)
  for (row_idx_reverse_iter l_ = z_nz.rbegin(); l_ != z_nz.rend(); l_++) {
    index_type l = *l_;
    for (row_iterator k_ = L.col(l).begin(); k_ != L.col(l).end(); k_++) {
      index_type k = p[ k_.index() ]; // indirection
      if (k < j) continue; // This is for the range (j:m)
      value_type L_kl = *k_;
      w[k] -= L_kl * w[l];
    }
  }
}

template<class Matrix>
inline
typename gilbert_peierls_lu_factorization<Matrix>::index_type
gilbert_peierls_lu_factorization<Matrix>::
find_pivot(index_type j)
{
  // Find the largest-magnitude element v(mu) of v(j:n)
  index_type mu = j;
  value_type largest = std::abs(w[j]);
  for (row_idx_reverse_iter l_ = v_nz.rbegin(); l_ != v_nz.rend(); l_++) {
    index_type l = *l_;
    value_type magnitude = std::abs(w[l]);
    if (magnitude > largest) {
      mu = l;
      largest = magnitude;
    }
  }
  return mu;
}

template<class Matrix>
inline
void gilbert_peierls_lu_factorization<Matrix>::
swap_rows(index_type i1, index_type i2) {
  // p := (i1,i2) o p
  p[ p_inv[i1] ]  = i2;
  p[ p_inv[i2]  ] = i1;

  // p^-1 := p^-1 o (i1,i2)
  std::swap(p_inv[i1], p_inv[i2]);

  std::swap(w[i1], w[i2]);
}

template<class Matrix>
inline
void gilbert_peierls_lu_factorization<Matrix>::
copy_z_into_U(index_type j)
{
  // construct the sparse U(0:j,j) from w and reset w to 0
  /* Only reset to 0 those elements which may be non-zero
  (same trick as for colour in the DFS) */
  for(row_idx_iter i_ = z_nz.begin(); i_ != z_nz.end(); i_++) {
    index_type i = *i_;
    U(i,j) = w[i];
    w[i] = 0;
  }
}

template<class Matrix>
inline
void gilbert_peierls_lu_factorization<Matrix>::
copy_v_into_L(index_type j)
{
  // construct the sparse L(:,j) from w and reset w to 0
  /* We need to store L unpermutated since we permute rows at each step.
  Same trick as copy_z_into_U to reset w. */
  double pivot = w[j];
  SCITBX_ASSERT(pivot);
  w[j] = 0;
  L(p_inv[j], j) = 1.;
  for(row_idx_reverse_iter i_ = v_nz.rbegin(); i_ != v_nz.rend(); i_++) {
    index_type i = *i_;
    if (i == j) continue;
    L(p_inv[i], j) = w[i]/pivot;
    w[i] = 0;
  }
}

}} // namespace scitbx::sparse

#endif
