#ifndef SCITBX_MATRIX_SYMMETRIC_RANK_1_UPDATE_H
#define SCITBX_MATRIX_SYMMETRIC_RANK_1_UPDATE_H

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/packed_matrix.h>
#include <scitbx/matrix/matrix_vector_operations.h>
#include <scitbx/matrix/vector_operations.h>
#include <scitbx/array_family/simple_io.h>
#include <fast_linalg/lapacke.h>
#include <algorithm>

namespace scitbx { namespace matrix {

  /// Sum of symmetric rank-1 updates \f$\alpha_i x_i x_i^T\f$
  template <typename T>
  class sum_of_symmetric_rank_1_updates
  {
  private:
    af::versa<T, af::packed_u_accessor> a;
    /// filled on first use by open_row(), and only if that is ever called
    af::shared<T> row_scratch;

  public:
    /// Initialise the sum to a zero matrix of size n
    /** The second argument is the row buffer size the other accumulator takes.
        This one keeps no buffer, so it is accepted and ignored, and the two
        remain interchangeable as template arguments.
     */
    sum_of_symmetric_rank_1_updates(int n, std::size_t=0)
    : a(n)
    {}

    /// Per-thread scratch the OpenMP accumulation needs for the matrix
    /** This one has no way to fold rows in by itself, so the OpenMP path does
        the packed rank-1 updates on its own, one matrix per thread.
     */
    static std::size_t omp_matrix_scratch(int n, int threads) {
      return std::size_t(threads)*(std::size_t(n)*(n + 1)/2);
    }

    /// Add \f$\alpha x x^T\f$ to the sum
    void add(af::const_ref<T> const &x, T alpha) {
      SCITBX_ASSERT(x.size() == a.accessor().n_rows())
      (x.size())
      (a.accessor().n_rows());
      add(x.begin(), alpha);
    }

    /// Add \f$\alpha x x^T\f$ (overload without size checks for speed)
    void add(T const *x, T alpha) {
      symmetric_packed_u_rank_1_update(a.accessor().n_rows(), a.begin(), x, alpha);
    }

    /** @brief A row for a caller which would otherwise build one and hand it
        over to add(), so that it may write it here directly instead.

    open_row() returns n_rows scalars to fill; commit_row() then does with them
    what add() would have done. The caller may read the row back between the two
    -- the least squares wants it for its own vector sums -- and it is untouched
    until commit_row.

    This one has nowhere to put a row but a scratch vector of its own, so a
    copy is still made; what is saved is the caller's copy, not this one. The
    buffered accumulator has somewhere, and saves both.
    */
    //@{
    T *open_row() {
      if (row_scratch.size() != a.accessor().n_rows()) {
        row_scratch.resize(a.accessor().n_rows());
      }
      return row_scratch.begin();
    }
    void commit_row(T alpha) { add(row_scratch.begin(), alpha); }
    //@}

    /// Add s in-place
    sum_of_symmetric_rank_1_updates
    &operator+=(sum_of_symmetric_rank_1_updates const &s) {
      af::ref<T> a_r = a.ref().as_1d();
      a_r += s.a.const_ref().as_1d();
      return *this;
    }

    /// Cancel all the rank-1 updates
    /** The sum is reset to the zero matrix */
    void reset() {
      std::fill(a.begin(), a.end(), T(0));
    }

    /// Called after after all rank-1 updates have been performed
    /** This does nothing in this class */
    void finalise() {}

    /// The resulting (symmetric) matrix
    /** It returns a meaningful result only after finalise() has been called */
    operator af::versa<T, af::packed_u_accessor>() {
      return a;
    }
  };

  /// Symmetric rank-N update \f$A^T A\f$, specified row by row
  /// but computed at BLAS 3 speed
  /** Thus this class is equivalent to class sum_of_symmetric_rank_1_updates
   *  when all \f$\alpha_i\f$ are non-negative, since then
   *  \f$\alpha x x^T = y y^T\f$ where \f$y = \sqrt{\alpha} x\f$ is one row
   *  of matrix A.
   *
   *  The rows are buffered and folded into the result a chunk at a time, so
   *  that the buffer is a working set rather than the whole of A. Holding all
   *  of A costs n_rows n_cols scalars, which on a large problem runs to several
   *  gigabytes and reallocates its way there, and there is no need to:
   *  \f$A^T A = \sum_c A_c^T A_c\f$ over any partition of the rows into chunks,
   *  and syrk accumulates into its output, so one call per chunk gives the same
   *  matrix. The flops are unchanged; what is paid for it is one pass over the
   *  result per chunk, which against the \f$n_{cols}^2 n_{rows}/2\f$ of the
   *  updates themselves is small.
   *
   *  Reassociating a floating-point sum is not bit-exact, so the matrix differs
   *  in the last digit or two from the unchunked one, as it does between any
   *  two summation orders.
   */
  template <typename T>
  class rank_n_update
  {
  public:
    /// Rows per chunk are chosen to fill about this much of a buffer
    /** A fallback only: a caller which knows its memory budget should say so,
     *  c.f. crystallographic_ls.accumulator_buffer_bytes. Generous, because
     *  each chunk costs a pass over the whole result and syrk's result is a
     *  full n^2 where the packed one it replaced was half that -- the buffer
     *  size hardly mattered with the packed update and does with this one, so
     *  the chunks are kept few. Still well below holding every row.
     *
     *  Functions rather than static constants so that taking their value never
     *  needs an out-of-line definition of them.
     */
    static std::size_t default_buffer_bytes() { return 256u << 20; }
    /// Below this, syrk is being called on too thin a matrix to pay for itself
    static std::size_t min_chunk_rows() { return 64; }

    /// CBLAS spells these as an enum which is not in scope here; the values are
    /// fixed by the CBLAS standard
    //@{
    static int cblas_row_major() { return 101; }
    static int cblas_upper() { return 121; }
    static int cblas_trans() { return 112; }
    //@}

    /// Prepare for a resulting matrix of size n
    /** buffer_bytes bounds the row buffer. A problem small enough that all of
     *  its rows fit is never chunked at all, so nothing about a small structure
     *  changes: it does the one update it always did.
     */
    rank_n_update(int n, std::size_t buffer_bytes=0)
    : aaT(std::size_t(n)*n, T(0)), aaT_packed(n), cols(n)
    {
      if (buffer_bytes == 0) {
        buffer_bytes = default_buffer_bytes();
      }
      std::size_t const row_bytes =
        std::max<std::size_t>(1, std::size_t(cols)*sizeof(T));
      chunk_rows = std::max(min_chunk_rows(), buffer_bytes/row_bytes);
      /* Not allocated here. The buffer grows geometrically to what is actually
         used and stops at one chunk, which is the point of chunking. Sizing it
         to the whole chunk up front measures markedly slower on a small problem
         -- where a generous budget makes a chunk very much larger than the whole
         of A -- with the row writes themselves proved free, so it is the buffer
         the update is handed that differs and not the filling of it.
       */
      reset();
    }

    /// Add a row \f$\sqrt{\alpha} x\f$ to matrix A
    /// Precondition: alpha >= 0
    void add(af::const_ref<T> const &x, T alpha) {
      SCITBX_ASSERT(x.size() == cols)(x.size())(cols);
      add(x.begin(), alpha);
    }

    /// Overload without size check for speed but alpha >= 0 still enforced
    /** Copy first, scale second, and not fused into one loop over the row.
     *  Fusing looks like it saves a pass and measures slower: a hand-written
     *  `dst[j] = x[j]*s` cannot be vectorised, the compiler having to allow for
     *  dst and x overlapping, where a copy and a scal are each free of that and
     *  each vectorised. The row is in L1 for the second pass anyway, having just
     *  been written, so there was nothing much to save.
     */
    void add(T const *x, T alpha) {
      T *row = open_row();
      std::copy(x, x + cols, row);
      commit_row(alpha);
    }

    /** @brief A row for a caller which would otherwise build one and hand it
        over to add(), so that it may write it here directly instead.

    This is the buffer the rows of A are already accumulated in, so a caller
    which fills the row itself removes a copy of it per equation: the gradient
    vector is built once, where before it was built and then copied in.

    The row comes back uninitialised -- the caller is expected to write all of
    it -- and stays readable and unscaled until commit_row, which the least
    squares depends on, its own vector sums reading the unweighted gradients.
    */
    //@{
    T *open_row() {
      if (n_pending_ == chunk_rows) {
        flush();
      }
      reserve_rows(n_pending_ + 1);
      return a.begin() + n_pending_*std::size_t(cols);
    }
    void commit_row(T alpha) {
      SCITBX_ASSERT(alpha >= 0)(alpha);
      matrix::scale_vector(cols, a.begin() + n_pending_*std::size_t(cols),
                           std::sqrt(alpha));
      n_pending_++;
    }
    //@}

    /// Add u in-place
    /** u is const, so its pending rows cannot be flushed through itself; they
     *  are folded straight into this one's result instead. Both that and the
     *  part u has already folded into its own are additions into the same
     *  matrix, so the order of the three below does not matter.
     */
    rank_n_update &operator+=(rank_n_update const &u) {
      SCITBX_ASSERT(u.cols == cols)(u.cols)(cols);
      flush();
      if (cols != 0 && u.n_pending_ != 0) {
        update(u.a.begin(), u.n_pending_);
      }
      if (u.folded) {
        // nothing has written this one's result yet, so u's is it rather than
        // something to add to whatever happens to be in there
        if (folded) {
          for (std::size_t i=0; i < aaT.size(); ++i) {
            aaT[i] += u.aaT[i];
          }
        }
        else {
          std::copy(u.aaT.begin(), u.aaT.end(), aaT.begin());
          folded = true;
        }
      }
      return *this;
    }

    /// Cancel all the rank-1 updates
    /** aaT_rfp is not cleared: the first chunk folded in overwrites it, which is
     *  both cheaper than zeroing it and what the unchunked version did.
     */
    void reset() {
      n_pending_ = 0;
      folded = false;
    }

    /** @brief The row buffer, for a caller which fills it itself.

    The OpenMP accumulation writes the scaled rows straight in, in parallel,
    rather than one at a time through add(). It needs no scratch matrix of its
    own: the syrk which folds a chunk in is threaded by the BLAS.
    */
    //@{
    static std::size_t omp_matrix_scratch(int, int) { return 0; }
    std::size_t n_pending() const { return n_pending_; }
    std::size_t max_pending() const { return chunk_rows; }
    /// Make room for n more rows and hand back the first of them
    T *open_pending(std::size_t n) {
      std::size_t const at = n_pending_;
      reserve_rows(at + n);
      n_pending_ += n;
      return a.begin() + at*std::size_t(cols);
    }
    void fold_pending() { flush(); }
    //@}

    /// Called after after all rank-1 updates have been performed
    void finalise() {
      if (cols == 0) {
        return;
      }
      flush();
      if (!folded) {
        // no rows at all: the result is the zero matrix, and nothing has
        // written the upper triangle for the copy below to read
        std::fill(aaT.begin(), aaT.end(), T(0));
        folded = true;
      }
      /* The upper triangle of a row-major square, read row by row, is the order
         packed_u stores it in, so this is one sequential pass and no
         arithmetic -- nothing here can lose a digit. The lower triangle syrk
         leaves untouched is never read; it is zero from construction so that
         operator+= may add the whole array without meeting whatever would
         otherwise be there.
       */
      T *p = aaT_packed.begin();
      for (int i=0; i < cols; ++i) {
        T const *row = &aaT[std::size_t(i)*cols];
        for (int j=i; j < cols; ++j) {
          *p++ = row[j];
        }
      }
    }

    /// The resulting (symmetric) matrix
    /** It returns a meaningful result only after finalise() has been called */
    operator af::versa<T, af::packed_u_accessor>() {
      return aaT_packed;
    }

  private:
    /// Fold the buffered rows into the result and empty the buffer
    void flush() {
      if (cols == 0 || n_pending_ == 0) {
        return;
      }
      update(a.begin(), n_pending_);
      n_pending_ = 0;
    }

    /// Make sure the buffer holds n rows, growing it geometrically if not
    /** Doubling, and never past one chunk, so that filling a chunk costs a
     *  logarithmic number of reallocations however the rows arrive. Growing it
     *  rather than sizing it in the constructor is what keeps a small problem
     *  from paying for a buffer sized by the memory budget instead of by its own
     *  number of reflections.
     */
    void reserve_rows(std::size_t n) {
      std::size_t const need = n*std::size_t(cols);
      if (a.size() >= need) {
        return;
      }
      std::size_t want = std::max(need, 2*a.size());
      want = std::min(want, chunk_rows*std::size_t(cols));
      a.resize(std::max(want, need));
    }

    /// aaT = A^T A for the first chunk, += it for the rest
    /** The buffer holds the rows of A, so row-major n_rows x cols, and
     *  C = A^T A is syrk with trans. Upper, because the upper triangle of a
     *  row-major square is laid out exactly as af::packed_u_accessor wants it,
     *  which makes finalise a straight copy.
     *
     *  syrk rather than the rectangular-full-packed sfrk, which was here
     *  before: sfrk halves the storage of the result and measures slower on
     *  the same shape, RFP being the less travelled path through the BLAS. The
     *  result costs n^2 rather than n(n+1)/2 for it.
     *
     *  beta = 1 accumulates the chunks; the first uses beta = 0 so that nothing
     *  need be zeroed beforehand.
     */
    void update(T const *rows, std::size_t n_rows) {
      using namespace fast_linalg;
      syrk(cblas_row_major(), cblas_upper(), cblas_trans(),
           cols, static_cast<int>(n_rows),
           1.0, rows, cols, folded ? 1.0 : 0.0, &aaT[0], cols);
      folded = true;
    }

    /// the rows of A not yet folded in; its size is capacity, n_pending_ is use
    af::shared<T> a;
    /// cols x cols, row major; syrk writes its upper triangle, the lower stays
    /// zero from construction
    af::shared<T> aaT;
    af::versa<T, af::packed_u_accessor> aaT_packed;
    int cols;
    std::size_t chunk_rows, n_pending_;
    /// whether anything has been folded into aaT_rfp since the last reset
    bool folded;
  };

}}


#endif
