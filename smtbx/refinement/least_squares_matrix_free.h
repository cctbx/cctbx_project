#ifndef SMTBX_REFINEMENT_LEAST_SQUARES_MATRIX_FREE_H
#define SMTBX_REFINEMENT_LEAST_SQUARES_MATRIX_FREE_H

/// Conjugate gradient least squares without ever forming the normal matrix.
/**
  The Gauss-Newton system of a refinement with a separable scale factor is
  A = D^T W D and b = D^T W r, where, with yc the computed observables, yo the
  observed, k* the scale factor the objective profiles out and all dots
  weighted,

      k*      = (w yo . yc)/(w yc . yc)
      r       = yo - k* yc
      grad_k* = (J^T W r - k* J^T W yc)/(w yc . yc)
      D       = k* J + yc (x) grad_k*

  Note it is D and not J: treating the system as J^T W J would be wrong by the
  rank-one term the scale factor contributes.

  Forming A costs O(n_par^2) storage, which is what stops a full-matrix
  refinement well before the structure factors become the difficulty. The
  smtbx::refinement::least_squares builders accumulate into whatever object is
  handed to them -- the accumulator is a template parameter of the constructor,
  not a fixed type -- so a class which accumulates something else entirely can
  reuse that machinery unchanged, reflection loop, twinning, solvent mask,
  Fc corrections, threading and all.

  Two such accumulators live here:

  - summary, which gathers in a single pass everything the conjugate gradients
    need to start: k*, grad k*, the right hand side, and the diagonal blocks
    that precondition them.

  - product, which computes A p in a single pass, and is what a conjugate
    gradient iteration calls.

  Both hold O(n_par) state per thread, against the O(n_par^2) a normal matrix
  needs, which is the entire point: the per-thread copies of the latter are
  what force the OpenMP path in least_squares_omp.h to chunk itself against
  max_memory: a packed normal matrix grows with the square of the parameter
  count before it is multiplied by the number of threads.
*/

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/ref.h>
#include <smtbx/error.h>

#include <boost/shared_ptr.hpp>
#include <vector>

namespace smtbx { namespace refinement { namespace least_squares {

  namespace af = scitbx::af;

  /// af::shared builds from a pair of pointers, which a vector's iterators
  /// are not obliged to be
  template <typename T>
  inline af::shared<T> as_shared(std::vector<T> const &v) {
    if (v.empty()) return af::shared<T>();
    return af::shared<T>(&v[0], &v[0] + v.size());
  }

  /// Which parameters are grouped together in the preconditioner.
  /** The blocks are the parameters of one scatterer -- its three coordinates
      and its six anisotropic displacement parameters -- taken together.
      Preconditioning with the diagonal alone, which is what Hendrickson and
      Konnert prescribe, leaves the condition number of these systems within a
      factor of two of doing nothing at all; blocking by scatterer drops it by
      two to three orders of magnitude and holds the conjugate gradients at
      fourteen to eighteen iterations whatever the size of the structure. It is
      not an optimisation, it is what makes the method work.

      The parameters of a block need not be contiguous, constraints being free
      to tie a scatterer's parameters to parameters elsewhere, so they are held
      as an explicit list of indices rather than as a range.
   */
  struct block_structure {
    /// the parameter indices of every block, one block after another
    std::vector<int> parameters;
    /// where each block starts in parameters, with a trailing total
    std::vector<int> parameter_start;
    /// where each block's dense n x n storage starts, with a trailing total
    std::vector<int> data_start;

    block_structure() {
      parameter_start.push_back(0);
      data_start.push_back(0);
    }

    block_structure(af::const_ref<int> const &parameters_,
                    af::const_ref<int> const &sizes)
    {
      parameters.assign(parameters_.begin(), parameters_.end());
      parameter_start.reserve(sizes.size() + 1);
      data_start.reserve(sizes.size() + 1);
      parameter_start.push_back(0);
      data_start.push_back(0);
      int at = 0, data_at = 0;
      for (std::size_t i = 0; i < sizes.size(); i++) {
        at += sizes[i];
        data_at += sizes[i]*sizes[i];
        parameter_start.push_back(at);
        data_start.push_back(data_at);
      }
      SMTBX_ASSERT(at == static_cast<int>(parameters.size()))
        (at)(parameters.size());
    }

    std::size_t n_blocks() const { return parameter_start.size() - 1; }
    int size_of(std::size_t i) const {
      return parameter_start[i + 1] - parameter_start[i];
    }
    int data_size() const { return data_start.back(); }
  };


  /// Common ground: the builders construct per-thread copies as T(n_parameters)
  /** The reflection loop hands each thread its own accumulator, made with that
      one constructor and merged afterwards with operator+=. Anything an
      accumulator needs beyond the number of parameters therefore cannot arrive
      through it, so it is left here for that constructor to pick up. A build is
      one call and the per-thread copies are all made on the calling thread
      before any of them runs, so this is only ever written and read from there.
   */
  template <typename Derived>
  struct thread_local_context {
    static boost::shared_ptr<Derived> &pending() {
      static boost::shared_ptr<Derived> instance;
      return instance;
    }

    /** @brief Per-thread scratch the OpenMP normal matrix accumulation needs.

    None: these accumulate a design matrix or a product with one and never form
    a normal matrix at all, so add_equations_omp is not a path they take. The
    builder asks whatever it was handed, so they have to be able to answer.
    */
    static std::size_t omp_matrix_scratch(int, int) { return 0; }
  };


  /// One pass gathering everything the conjugate gradients need to start.
  template <typename FloatType>
  class separable_scale_factor_summary
    : public thread_local_context<separable_scale_factor_summary<FloatType> >
  {
  public:
    typedef FloatType scalar_t;
    typedef thread_local_context<separable_scale_factor_summary<FloatType> >
      context_t;

    /// The one Python builds; it leaves itself where the per-thread copies look
    separable_scale_factor_summary(int n_parameters,
                                   af::const_ref<int> const &block_parameters,
                                   af::const_ref<int> const &block_sizes)
      : n_params(n_parameters),
        blocks_(new block_structure(block_parameters, block_sizes))
    {
      allocate();
      context_t::pending().reset(new separable_scale_factor_summary(*this));
    }

    /// The one the reflection loop makes for each thread
    separable_scale_factor_summary(int n_parameters)
      : n_params(n_parameters)
    {
      boost::shared_ptr<separable_scale_factor_summary> const &source
        = context_t::pending();
      SMTBX_ASSERT(source);
      SMTBX_ASSERT(source->n_params == n_parameters)
        (source->n_params)(n_parameters);
      blocks_ = source->blocks_;
      allocate();
    }

    int n_parameters() const { return n_params; }

    void add_residual(scalar_t yc, scalar_t yo, scalar_t w) {
      sum_w_yo_sq += w*yo*yo;
      sum_w_yc_sq += w*yc*yc;
      sum_w_yo_yc += w*yo*yc;
      n_equations++;
    }

    void add_equation(scalar_t yc, af::const_ref<scalar_t> const &grad_yc,
                      scalar_t yo, scalar_t w)
    {
      add_residual(yc, yo, w);
      scalar_t const w_yo = w*yo, w_yc = w*yc;
      for (int i = 0; i < n_params; i++) {
        yo_dot_grad_yc[i] += w_yo*grad_yc[i];
        yc_dot_grad_yc[i] += w_yc*grad_yc[i];
      }
      // the diagonal blocks of J^T W J; the scale factor turns them into the
      // blocks of D^T W D in finalise, once k* and grad k* are known
      block_structure const &b = *blocks_;
      for (std::size_t bi = 0; bi < b.n_blocks(); bi++) {
        int const n = b.size_of(bi);
        int const *p = &b.parameters[b.parameter_start[bi]];
        gathered.resize(n);
        for (int a = 0; a < n; a++) gathered[a] = grad_yc[p[a]];
        scalar_t *data = &block_data[b.data_start[bi]];
        for (int a = 0; a < n; a++) {
          scalar_t const row = w*gathered[a];
          for (int c = 0; c < n; c++) data[a*n + c] += row*gathered[c];
        }
      }
    }

    separable_scale_factor_summary &
    operator+=(separable_scale_factor_summary const &other) {
      SMTBX_ASSERT(n_params == other.n_params)(n_params)(other.n_params);
      sum_w_yo_sq += other.sum_w_yo_sq;
      sum_w_yc_sq += other.sum_w_yc_sq;
      sum_w_yo_yc += other.sum_w_yo_yc;
      n_equations += other.n_equations;
      for (int i = 0; i < n_params; i++) {
        yo_dot_grad_yc[i] += other.yo_dot_grad_yc[i];
        yc_dot_grad_yc[i] += other.yc_dot_grad_yc[i];
      }
      for (std::size_t i = 0; i < block_data.size(); i++) {
        block_data[i] += other.block_data[i];
      }
      return *this;
    }

    /// Assemble what was accumulated into what the solver asks for.
    void finalise(bool objective_only = false) {
      SMTBX_ASSERT(sum_w_yc_sq > 0)(sum_w_yc_sq);
      scale_factor_ = sum_w_yo_yc/sum_w_yc_sq;
      scalar_t const k = scale_factor_;
      // sum w r^2 = sum w yo^2 - k sum w yo yc, the cross terms cancelling
      // because k is the value which minimises it
      objective_ = (sum_w_yo_sq - k*sum_w_yo_yc)/(2*sum_w_yo_sq);
      if (objective_only) {
        finalised = true;
        return;
      }
      grad_scale_factor_.resize(n_params);
      right_hand_side_.resize(n_params);
      for (int i = 0; i < n_params; i++) {
        // grad k* = (J^T W r - k J^T W yc)/(w yc . yc), with J^T W r expanded
        // as J^T W yo - k J^T W yc
        grad_scale_factor_[i] =
          (yo_dot_grad_yc[i] - 2*k*yc_dot_grad_yc[i])/sum_w_yc_sq;
        // b = D^T W r = k J^T W r + (w r . yc) grad k*, and the second term is
        // identically zero: w r . yc = sum w yo yc - k sum w yc^2 = 0
        right_hand_side_[i] =
          k*(yo_dot_grad_yc[i] - k*yc_dot_grad_yc[i])/sum_w_yo_sq;
      }
      // D[b] (x) D[b] = k^2 g (x) g + k (yc.g (x) grad k + grad k (x) yc.g)
      //                 + (w yc . yc) grad k (x) grad k
      block_structure const &b = *blocks_;
      for (std::size_t bi = 0; bi < b.n_blocks(); bi++) {
        int const n = b.size_of(bi);
        int const *p = &b.parameters[b.parameter_start[bi]];
        scalar_t *data = &block_data[b.data_start[bi]];
        for (int a = 0; a < n; a++) {
          scalar_t const yc_a = yc_dot_grad_yc[p[a]], gk_a =
            grad_scale_factor_[p[a]];
          for (int c = 0; c < n; c++) {
            scalar_t const yc_c = yc_dot_grad_yc[p[c]], gk_c =
              grad_scale_factor_[p[c]];
            data[a*n + c] = (k*k*data[a*n + c] + k*(yc_a*gk_c + gk_a*yc_c)
                             + sum_w_yc_sq*gk_a*gk_c)/sum_w_yo_sq;
          }
        }
      }
      finalised = true;
    }

    scalar_t scale_factor() const { check(); return scale_factor_; }
    scalar_t objective() const { check(); return objective_; }
    scalar_t sum_w_yo_sq_value() const { return sum_w_yo_sq; }
    scalar_t sum_w_yc_sq_value() const { return sum_w_yc_sq; }
    int n_equations_value() const { return n_equations; }

    af::shared<scalar_t> grad_scale_factor() const {
      check();
      return as_shared(grad_scale_factor_);
    }
    af::shared<scalar_t> right_hand_side() const {
      check();
      return as_shared(right_hand_side_);
    }
    /// the blocks, one dense n x n after another in the order they were given
    af::shared<scalar_t> blocks() const {
      check();
      return as_shared(block_data);
    }

    /// the OpenMP path allocates a normal matrix per thread before it consults
    /// the accumulator at all, so it is no use here; see the header comment
    template <class T1, class T2, class T3>
    void add_residuals_omp(int, int, int, T1 const &, T2 const &, T3 const &) {
      SMTBX_NOT_IMPLEMENTED();
    }
    template <class T1, class T2>
    void add_residuals_omp(int, int, int, T1 const &, T2 const &) {
      SMTBX_NOT_IMPLEMENTED();
    }
    template <class T1, class T2, class T3, class T4, class T5, class T6,
              class T7>
    void add_equations_omp(int, int, int, int, int, T1 &, T2 &, T3 &, T4 const &,
                           T5 &, T6 const &, T7 const &) {
      SMTBX_NOT_IMPLEMENTED();
    }

  private:
    void allocate() {
      yo_dot_grad_yc.assign(n_params, 0);
      yc_dot_grad_yc.assign(n_params, 0);
      block_data.assign(blocks_->data_size(), 0);
      sum_w_yo_sq = sum_w_yc_sq = sum_w_yo_yc = 0;
      scale_factor_ = objective_ = 0;
      n_equations = 0;
      finalised = false;
    }
    void check() const { SMTBX_ASSERT(finalised); }

    int n_params;
    boost::shared_ptr<block_structure> blocks_;
    std::vector<scalar_t> yo_dot_grad_yc, yc_dot_grad_yc, block_data;
    std::vector<scalar_t> grad_scale_factor_, right_hand_side_;
    mutable std::vector<scalar_t> gathered;
    scalar_t sum_w_yo_sq, sum_w_yc_sq, sum_w_yo_yc;
    scalar_t scale_factor_, objective_;
    int n_equations;
    bool finalised;
  };


  /// One pass computing A p, with the linearisation held fixed.
  /** The structure does not move during the conjugate gradients, so yc, the
      weights, k* and grad k* are constants across the iterations of a cycle and
      are supplied rather than recomputed. What has to be recomputed is the
      gradient of each reflection, which is the expensive part and the price of
      not storing it.

          for each reflection h:
              g   = grad yc_h
              t   = w_h (k (g . p) + yc_h (grad k . p))    ( = w_h D_h . p )
              y  += t g
              s  += t yc_h
          A p = (k y + s grad k)/sum w yo^2
   */
  template <typename FloatType>
  class separable_scale_factor_product
    : public thread_local_context<separable_scale_factor_product<FloatType> >
  {
  public:
    typedef FloatType scalar_t;
    typedef thread_local_context<separable_scale_factor_product<FloatType> >
      context_t;

    separable_scale_factor_product(int n_parameters,
                                   scalar_t scale_factor,
                                   af::const_ref<scalar_t> const &grad_k,
                                   af::const_ref<scalar_t> const &p,
                                   scalar_t sum_w_yo_sq)
      : n_params(n_parameters),
        k(scale_factor),
        sum_w_yo_sq_(sum_w_yo_sq),
        grad_k_(grad_k.begin(), grad_k.end()),
        p_(p.begin(), p.end())
    {
      SMTBX_ASSERT(grad_k.size() == n_parameters)(grad_k.size())(n_parameters);
      SMTBX_ASSERT(p.size() == n_parameters)(p.size())(n_parameters);
      grad_k_dot_p = 0;
      for (int i = 0; i < n_params; i++) grad_k_dot_p += grad_k_[i]*p_[i];
      allocate();
      context_t::pending().reset(new separable_scale_factor_product(*this));
    }

    separable_scale_factor_product(int n_parameters)
      : n_params(n_parameters)
    {
      boost::shared_ptr<separable_scale_factor_product> const &source
        = context_t::pending();
      SMTBX_ASSERT(source);
      SMTBX_ASSERT(source->n_params == n_parameters)
        (source->n_params)(n_parameters);
      k = source->k;
      sum_w_yo_sq_ = source->sum_w_yo_sq_;
      grad_k_ = source->grad_k_;
      p_ = source->p_;
      grad_k_dot_p = source->grad_k_dot_p;
      allocate();
    }

    int n_parameters() const { return n_params; }

    void add_residual(scalar_t, scalar_t, scalar_t) {}

    void add_equation(scalar_t yc, af::const_ref<scalar_t> const &grad_yc,
                      scalar_t, scalar_t w)
    {
      scalar_t g_dot_p = 0;
      for (int i = 0; i < n_params; i++) g_dot_p += grad_yc[i]*p_[i];
      scalar_t const t = w*(k*g_dot_p + yc*grad_k_dot_p);
      for (int i = 0; i < n_params; i++) y_[i] += t*grad_yc[i];
      s_ += t*yc;
    }

    separable_scale_factor_product &
    operator+=(separable_scale_factor_product const &other) {
      SMTBX_ASSERT(n_params == other.n_params)(n_params)(other.n_params);
      for (int i = 0; i < n_params; i++) y_[i] += other.y_[i];
      s_ += other.s_;
      return *this;
    }

    void finalise(bool objective_only = false) {
      if (objective_only) return;
      for (int i = 0; i < n_params; i++) {
        y_[i] = (k*y_[i] + s_*grad_k_[i])/sum_w_yo_sq_;
      }
      finalised = true;
    }

    af::shared<scalar_t> result() const {
      SMTBX_ASSERT(finalised);
      return as_shared(y_);
    }

    template <class T1, class T2, class T3>
    void add_residuals_omp(int, int, int, T1 const &, T2 const &, T3 const &) {
      SMTBX_NOT_IMPLEMENTED();
    }
    template <class T1, class T2>
    void add_residuals_omp(int, int, int, T1 const &, T2 const &) {
      SMTBX_NOT_IMPLEMENTED();
    }
    template <class T1, class T2, class T3, class T4, class T5, class T6,
              class T7>
    void add_equations_omp(int, int, int, int, int, T1 &, T2 &, T3 &, T4 const &,
                           T5 &, T6 const &, T7 const &) {
      SMTBX_NOT_IMPLEMENTED();
    }

  private:
    void allocate() {
      y_.assign(n_params, 0);
      s_ = 0;
      finalised = false;
    }

    int n_params;
    scalar_t k, sum_w_yo_sq_, grad_k_dot_p, s_;
    std::vector<scalar_t> grad_k_, p_, y_;
    bool finalised;
  };

}}}

#endif // GUARD
