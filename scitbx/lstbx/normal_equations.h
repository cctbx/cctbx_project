/// Tools implementing the Gauss-Newton method for non-linear least-squares.

#ifndef SCITBX_GAUSS_NEWTON_H
#define SCITBX_GAUSS_NEWTON_H

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/shared_algebra.h>
#include <scitbx/array_family/ref_algebra.h>
#include <scitbx/array_family/owning_ref.h>
#include <scitbx/array_family/accessors/row_and_column.h>
#include <scitbx/matrix/cholesky.h>
#include <scitbx/matrix/symmetric_rank_1_update.h>
#include <scitbx/sparse/matrix.h>
#include <scitbx/sparse/triangular.h>
#include <sstream>
#include <cmath>
#include <vector>
#if defined(_OPENMP)
  #include <omp.h>
#endif


namespace scitbx { namespace lstbx { namespace normal_equations {

#define SCITBX_LSTBX_DECLARE_ARRAY_TYPE(FloatType)                            \
    typedef FloatType scalar_t;                                               \
    typedef af::ref_owning_versa<scalar_t,                                    \
                                 af::packed_u_accessor>                       \
            symmetric_matrix_owning_ref_t;                                    \
    typedef af::ref_owning_versa<scalar_t,                                    \
                                 af::packed_u_accessor>                       \
            upper_diagonal_matrix_owning_ref_t;                               \
    typedef af::ref<scalar_t,                                                 \
                    af::packed_u_accessor>                                    \
            symmetric_matrix_ref_t;                                           \
    typedef af::versa<scalar_t,                                               \
                      af::packed_u_accessor>                                  \
            symmetric_matrix_t;                                               \
    typedef af::versa<scalar_t,                                               \
                      af::packed_u_accessor>                                  \
            upper_diagonal_matrix_t;                                          \
    typedef af::ref_owning_versa<FloatType, af::mat_grid> matrix_owning_ref_t;\
    typedef af::ref<FloatType, af::mat_grid> matrix_ref_t;                    \
    typedef af::ref_owning_shared<scalar_t> vector_owning_ref_t;              \
    typedef af::shared<scalar_t> vector_t;                                    \
    typedef af::ref<scalar_t> vector_ref_t;


  /// Normal equations for linear least-squares problem.
  /** The least-squares target reads

      \f[ L(x) = \| A x - b \|^2 \f]

      where the norm is diagonal-weighted

      \f[ \| y \|^2 = \sum_i w_i y_i^2 \f]

      Objects of this type may also be used to hold the normal equations
      from a non-linear problem after they have been built.
  */
  template <typename FloatType>
  class linear_ls
  {
  public:
    SCITBX_LSTBX_DECLARE_ARRAY_TYPE(FloatType);

    /// Construct a least-squares problem with the given number of unknowns.
    linear_ls(int n_parameters)
      : solved_(false),
        n_accumulated_(0),
        normal_matrix_(n_parameters),
        right_hand_side_(n_parameters)
    {}

    /// Number of unknown parameters
    int n_parameters() const { return right_hand_side_.size(); }

    /// Initialise the least-squares problem with the given normal matrix A
    /// and right hand side b
    linear_ls(symmetric_matrix_t const &a, vector_t const &b)
      : solved_(false),
        // Handed a matrix outright: it counts as accumulated, or the guard
        // below would treat a caller-supplied system as "never built".
        n_accumulated_(1),
        normal_matrix_(a),
        right_hand_side_(b)
    {
      SCITBX_ASSERT(a.accessor().n == b.size());
    }

    /// Add the equation \f$ A_{i.} x = b_i \f$ with the given weight
    void add_equation(scalar_t b_i,
                      af::const_ref<scalar_t> const &a_row,
                      scalar_t w)
    {
      ++n_accumulated_;
      scalar_t *p = normal_matrix_.begin();
      for (int i=0; i<n_parameters(); ++i)  {
        right_hand_side_[i] += w * a_row[i] * b_i;
        for (int j=i; j<n_parameters(); ++j) *p++ += w * a_row[i] * a_row[j];
      }
    }

    /// Add the equations A x = b with the given weights
    /** w[i] weights the i-th equation, i.e. the row \f$ A_{i.} \f$.
        If negate_right_hand_side, then the equation is A x + b = 0 instead
        Optimise_for_sparse can be used to control the method used to calculate
        A^T W A, which may vary significantly in performance depending on the
        shape and sparsity of the problem. Simple testing suggests that this
        variable should be set to true for a highly sparse case, but it is
        recommended that a developer assesses the best option for their problem.
        See github pull request #295 for further discussion on this topic.
     */
    void add_equations(af::const_ref<scalar_t> const &b,
                       sparse::matrix<scalar_t> const &a,
                       af::const_ref<scalar_t> const &w,
                       bool negate_right_hand_side=false,
                       bool optimise_for_tall_matrix=true)
    {
      SCITBX_ASSERT(   a.n_rows() == b.size()
                    && b.size()   == w.size())(a.n_rows())(b.size())(w.size());
      SCITBX_ASSERT(a.n_cols() == n_parameters());
      sparse::matrix<scalar_t> at_w_a;
      if (optimise_for_tall_matrix) {
        at_w_a = a.this_transpose_times_diagonal_times_this(w);
      } else {
        at_w_a = a.transpose().this_times_diagonal_times_this_transpose(w);
      }
      vector_t a_t_w_b = a.transpose_times((w * b).const_ref());
      update_matrix(at_w_a, a_t_w_b, negate_right_hand_side);
    }

    // Add directly to the normal matrix equation
    void update_matrix(sparse::matrix<scalar_t> const &at_w_a,
                       vector_t const &a_t_w_b,
                       bool negate_right_hand_side){
      ++n_accumulated_;
      normal_matrix_ += sparse::upper_diagonal_of(at_w_a);
      if (negate_right_hand_side) right_hand_side_ -= a_t_w_b.const_ref();
      else                        right_hand_side_ += a_t_w_b.const_ref();

    }

    /// Reset the state to construction time, i.e. no equations accumulated
    void reset() {
      solved_ = false;
      n_accumulated_ = 0;
      std::fill(normal_matrix_.begin(), normal_matrix_.end(), scalar_t(0));
      std::fill(right_hand_side_.begin(), right_hand_side_.end(), scalar_t(0));
    }

    /// Only available if the equations have not been solved yet
    symmetric_matrix_t normal_matrix() const {
      SCITBX_ASSERT(!solved());
      return normal_matrix_.array();
    }

    /// Only available if the equations have not been solved yet
    vector_t right_hand_side() const {
      SCITBX_ASSERT(!solved());
      return right_hand_side_.array();
    }

    /** \brief Solve the normal equations for the parameters (linear case)
         or their shift (linearised non-linear case)
     */
    void solve() {
      using scitbx::matrix::cholesky::u_transpose_u_decomposition_in_place;
      int const n = n_parameters();

      /* Jacobi (symmetric diagonal) preconditioning, i.e. decompose
         \f$ A' = D^{-1} A D^{-1} \f$ with \f$ D = \mathrm{diag}(\sqrt{A_{ii}}) \f$
         rather than A itself.

         Why: crystallographic normal matrices routinely span many orders of
         magnitude between parameter types, and the decomposition's backward
         error goes like \f$ n \epsilon \lambda_{max} \f$ -- a property of the
         *largest* eigenvalue. When the weakest parameter's eigenvalue is
         smaller than that, whether a pivot stays positive is decided by
         rounding rather than by the data, so the same model succeeds on one
         build and throws "not positive definite" on another.

         Measured on an electron-diffraction N-beam refinement (131 parameters,
         a refined sample thickness among them): lambda_min 8.5e-08 against a
         rounding floor of 7.5e-08 -- a margin of 1.1. Scaling takes the
         effective condition number from 3.0e13 to about 9, and the margin with
         it. The thickness row is 1e-07 where the ADP rows are 1e+06; that
         spread is a choice of units, not a statement about the data, and the
         decomposition should not be sensitive to it.

         This is *not* damping: A' and A have the same solution in exact
         arithmetic, no term is added, and nothing is regularised away. A
         genuinely undetermined parameter still fails, at the same index, which
         is what makes the error message worth reading.
       */
      /* Parameters no observation touches are frozen rather than fought over.

         For a normal matrix \f$ A = M^T W M \f$ with positive weights,
         \f$ A_{ii} = \sum_k w_k M_{ki}^2 = 0 \f$ forces column i of the design
         matrix to vanish identically -- so the whole row, the whole column and
         the right-hand side entry are zero too. Such a parameter is not merely
         ill-determined, it is absent from the system: freezing it at a zero
         shift provably cannot change any other parameter's shift, because
         nothing is coupled to it.

         Without this the decomposition fails at that pivot and the whole cycle
         is lost, reporting only an index. That has been a recurring failure in
         both the X-ray and the electron-diffraction paths, because there are
         many ways to end up with a live parameter nothing constrains: a
         rotatable AFIX group whose hydrogens are all fixed, an atom at zero
         occupancy, a parameter left refining after everything it acts on has
         been fixed.

         `zero_diagonal_indices()` reports what was frozen so a caller can say
         so. Silence would be worse than the crash: a user is entitled to know
         that a parameter they asked to refine carried no information.
       */
      zero_diagonal_indices_.clear();
      std::vector<scalar_t> d(n, scalar_t(1));
      {
        scalar_t const *p = normal_matrix_.begin();
        for (int i=0; i<n; ++i) {
          scalar_t a_ii = *p;
          if (a_ii > 0) d[i] = std::sqrt(a_ii);
          else zero_diagonal_indices_.push_back(i);
          p += (n - i);            // start of the next packed row
        }
      }
      /* The indices are *recorded*, not acted on: an unconstrained parameter
         still fails the decomposition, exactly as before.

         Freezing them at a zero shift was tried and reverted. It is sound in
         isolation -- for \f$ A = M^T W M \f$ a zero diagonal forces the whole
         row, column and rhs to vanish, so the parameter is decoupled -- but it
         cannot be told apart from the state CGLS leaves behind when it
         declines to build normal equations. There, restraints and origin
         fixing still accumulate (`cgls.py` calls
         `build_up(objective_only=True)`), so "was anything accumulated" does
         not separate the two either; freezing then manufactures an identity
         and hands back a covariance matrix describing no model, which
         `tst_cgls.exercise_without_standard_uncertainties` rightly rejects.

         It would also not have helped the case that prompted it. Olex2's
         electron-diffraction dump for a failing thickness step shows a 1x1
         matrix with diagonal 0 **and gradient 0**: the derivative never
         reached the normal equations at all. Letting that "succeed" would
         refine nothing while reporting success, which is worse than the
         failure it replaces.
       */

      {
        scalar_t *p = normal_matrix_.begin();
        for (int i=0; i<n; ++i) {
          right_hand_side_[i] /= d[i];
          for (int j=i; j<n; ++j) *p++ /= (d[i]*d[j]);
        }
      }

      u_transpose_u_decomposition_in_place<scalar_t> cholesky(normal_matrix_);
      if(cholesky.failure) {
        std::ostringstream buffer;
        buffer << "SCITBX_ASSERT(!cholesky.failure) failure in index: "
          << cholesky.failure.index;
        throw SCITBX_ERROR(buffer.str());
      }
      SCITBX_ASSERT(!cholesky.failure);
      cholesky.solve_in_place(right_hand_side_);

      /* Undo the scaling, in both things a caller can still reach.

         The solution first: solving A' y = b' gives \f$ y = D x \f$, so the
         step is \f$ x = D^{-1} y \f$.

         Then the factor. `cholesky_factor()` stays live after solve() and
         `smtbx.refinement.least_squares.covariance_matrix()` inverts it for
         every ESD the user ever sees, so it must be the factor of A and not of
         A'. From \f$ U'^T U' = D^{-1} A D^{-1} \f$ it follows that
         \f$ U = U' D \f$ satisfies \f$ U^T U = A \f$ -- one multiplication per
         stored element, by the scale of its *column*. Without this the ESDs
         would be wrong by \f$ d_i d_j \f$ and nothing would announce it.
       */
      for (int i=0; i<n; ++i) right_hand_side_[i] /= d[i];
      {
        scalar_t *p = normal_matrix_.begin();
        for (int i=0; i<n; ++i) {
          for (int j=i; j<n; ++j) *p++ *= d[j];
        }
      }
      solved_ = true;
    }

    bool solved() const {
      return solved_;
    }

    /// Only available after the equations have been solved
    upper_diagonal_matrix_t cholesky_factor() const {
      SCITBX_ASSERT(solved());
      return normal_matrix_.array();
    }

    /// Only available after the equations have been solved
    vector_t solution() const {
      SCITBX_ASSERT(solved());
      return right_hand_side_.array();
    }

    /// Parameters frozen by the last solve() because nothing constrained them
    /** Empty in a healthy refinement. A non-empty result names parameters the
        user asked to refine and no observation touched -- worth reporting.
     */
    af::shared<int> zero_diagonal_indices() const {
      af::shared<int> result;
      for (std::size_t i=0; i<zero_diagonal_indices_.size(); ++i) {
        result.push_back(zero_diagonal_indices_[i]);
      }
      return result;
    }

  public:
    bool solved_;
    std::size_t n_accumulated_;
    symmetric_matrix_owning_ref_t normal_matrix_;
    vector_owning_ref_t right_hand_side_;
    std::vector<int> zero_diagonal_indices_;
  };


  /// Normal equations for non-linear least-squares
  /** The least-squares target reads

      \f[ L(x) = \frac{1}{2} \|r(x)\|^2 \f]

      where the norm is diagonal-weighted

      \f[ \| y \|^2 = \sum_i w_i y_i^2 \f]

      and \f$r(x)\f$ is a vector of residuals depending on a vector
      of unknowns \f$x\f$.
  */
  template <typename FloatType>
  class non_linear_ls
  {
  public:
    SCITBX_LSTBX_DECLARE_ARRAY_TYPE(FloatType);

    /// Construct a least-squares problem with the given number of unknowns.
    non_linear_ls(int n_parameters)
      : n_equations_(0),
        r_sq(0),
        linearised(n_parameters)
    {}

    /// Construct with an exiting L.S. problem.
    /** That is
          - this->objective()  == objective
          - this->step_equations().right_hand_side() == opposite_of_grad_objective
          - this->step_equations().normal_matrix() == normal_matrix
     */
    non_linear_ls(std::size_t n_equations,
                  scalar_t objective,
                  vector_t const &opposite_of_grad_objective,
                  symmetric_matrix_t const &normal_matrix)
    : n_equations_(n_equations),
      r_sq(2*objective),
      linearised(normal_matrix, opposite_of_grad_objective)
    {}

    /// Number of equations
    /** i.e. number of components of the residual vector \f$r(x)\f$
     */
    std::size_t n_equations() const { return n_equations_; }

    /// Number of unknown parameters
    int n_parameters() const { return linearised.n_parameters(); }

    /// Number of degrees of freedom
    std::size_t dof() const { return n_equations() - n_parameters(); }

    /// Add the given residual with the given weight
    void add_residual(scalar_t r, scalar_t w) {
      n_equations_++;
      r_sq += w*r*r;
    }

    /// Add the given residuals with the given weights
    void add_residuals(af::const_ref<scalar_t> const &r,
                       af::const_ref<scalar_t> const &w)
    {
      for (int i=0; i<r.size(); ++i) {
        add_residual(r[i], w.size() ? w[i] : 1);
      }
    }

    /// Add the linearisation of the equation \f$r_i(x) = 0\f$
    /** with the given weight */
    void add_equation(scalar_t r,
                      af::const_ref<scalar_t> const &grad_r,
                      scalar_t w)
    {
      add_residual(r, w);
      linearised.add_equation(-r, grad_r, w);
    }

    /// Add the linearisation of the equations \f$r(x) = 0\f$ all at once
    /** The Jacobian is that of \f$x \mapto r(x)\f$.
     */
    void add_equations(af::const_ref<scalar_t> const &r,
                       af::const_ref<scalar_t, af::mat_grid> const &jacobian,
                       af::const_ref<scalar_t> const &w)
    {
      SCITBX_ASSERT(   r.size() == jacobian.n_rows()
                    && (!w.size() || r.size() == w.size()))
                   (r.size())(jacobian.n_rows())(w.size());
      SCITBX_ASSERT(jacobian.n_columns() == n_parameters())
                   (jacobian.n_columns())(n_parameters());
      for (int i=0; i<r.size(); ++i) {
        add_equation(r[i], af::row(jacobian, i), w.size() ? w[i] : 1);
      }
    }

    void add_equations(af::const_ref<scalar_t> const &r,
                       sparse::matrix<scalar_t> const &jacobian,
                       af::const_ref<scalar_t> const &w,
                       bool negate_right_hand_side=true,
                       bool optimise_for_tall_matrix=true)
    {
      SCITBX_ASSERT(   r.size() == jacobian.n_rows()
                    && (!w.size() || r.size() == w.size()))
                   (r.size())(jacobian.n_rows())(w.size());
      SCITBX_ASSERT(jacobian.n_cols() == n_parameters())
                   (jacobian.n_cols())(n_parameters());
      add_residuals(r, w);
      linearised.add_equations(r, jacobian, w, negate_right_hand_side, optimise_for_tall_matrix);
    }

    /// Objective value \f$L(x)\f$ for the current value of the unknowns
    scalar_t objective() const { return r_sq/2; }

    /// The \f$chi^2\f$ of the fit
    /**
        \f [ \frac{\sum_i w_i r_i(x)^2}
                  {n_{\text{equations}} - n_{\text{parameters}}

        Strictly speaking, this is only meaningful when the residuals have
        the form used in a fit, \f$r_i(x) = \text{model} - \text{data}\f$,
        but the computation is the same in the general case.
     */
    scalar_t chi_sq() const { return r_sq/dof(); }

    /// Linearised equations to solve for a step
    linear_ls<scalar_t> &step_equations() { return linearised; }

    /// Reset the state to construction time, i.e. no equations accumulated
    void reset() {
      n_equations_ = 0;
      r_sq = 0;
      linearised.reset();
    }

  protected:
    std::size_t n_equations_;
    scalar_t r_sq;
    linear_ls<scalar_t> linearised;
  };


  /// Normal equations for least-squares fit with an overall scale.
  /** The least-squares target reads

      \f[ L(K, x) = \frac{1}{2} \frac{ \sum w ( K y_c(x) - y_o )^2 }
                                     { \sum w y_o^2 }
      \f]

      where the both of \f$ y_c(x) \f$ and \f$ y_o \f$ are vectors,
      respectively the model to fit to the data. Alternatively, the
      non-normalised

      \f[ \tilde{L}(K, x) = \frac{1}{2} \sum w ( K y_c(x) - y_o )^2 \f]

      may be used instead.

      One takes advantage of the separability of the problem:

        - step 1: \f$ K^*(x) = \argmin_K L(K, x) \f$

        - step 2: Build the Newton equations for the problem
                  \f$ \min_x L(K^*(x), x) \f$

          in the Gauss approximation of small residuals (reduced equations).

  Reference:
   Separable nonlinear least squares
   H.B. Nielsen
   Technical report IMM-REP-2000-01
   http:http://www2.imm.dtu.dk/pubdb/views/edoc_download.php/646/ps/imm646.ps

   and references therein.
  */
  template <typename FloatType,
            template<typename> class SumOfRank1Updates=matrix::sum_of_symmetric_rank_1_updates>
  class non_linear_ls_with_separable_scale_factor
  {
  public:
    SCITBX_LSTBX_DECLARE_ARRAY_TYPE(FloatType);
    typedef SumOfRank1Updates<FloatType> sum_of_rank_1_updates_t;

    /// Construct a least-squares problem with the given number of parameters.
    /** That is the length of the vector \f$ x \f$. The flag normalised
        specify whether to use the normalised objective \f$L\f$ or the
        non-normalised objective \f$\tilde{L}\f$.
     */
    /** accumulator_buffer_bytes is passed on to the accumulator, for those of
        them which buffer rows before folding them into the normal matrix; zero
        leaves the choice to it. It has no meaning for an accumulator which does
        not buffer, and is ignored there.
     */
    non_linear_ls_with_separable_scale_factor(
      int n_parameters,
      bool normalised=true,
      std::size_t accumulator_buffer_bytes=0)
      : yo_dot_yc(0), yc_sq(0), yo_sq(0),
        n_params(n_parameters),
        n_data(0),
        normalised_(normalised),
        grad_yc_dot_grad_yc(n_parameters, accumulator_buffer_bytes),
        yo_dot_grad_yc(n_parameters),
        yc_dot_grad_yc(n_parameters),
        grad_k_star(n_parameters),
        finalised_(false),
        reduced_ls(n_parameters)
    {}

    /// Number of unknown parameters, not including the overall scale factor
    int n_parameters() const { return n_params; }

    /// Number of equations \f$y_o = K y_c(x)\f$ plus those added to
    /// the reduced_problem().
    std::size_t n_equations() const {
      return finalised() ? reduced_ls.n_equations() : n_data;
    }

    /// Number of degrees of freedom.
    /** This does take into account the equations added to the reduced_problem().
     */
    std::size_t dof() const { return n_equations() - n_parameters(); }

    /// Whether the L.S. target is normalised by \f$ \sum w y_o^2 \f$ or not
    bool normalised() const { return normalised_; }

    void add_residual(scalar_t yc, scalar_t yo, scalar_t w) {
      n_data++;
      yo_sq += w * yo * yo;
      yo_dot_yc += w * yo * yc;
      yc_sq += w * yc * yc;
    }

    /** \brief Add the linearisation of the equation
         \f$y_{c,i} \propto y_{o,i}\f$ with weight w.
     */
    void add_equation(scalar_t yc, af::const_ref<scalar_t> const &grad_yc,
                      scalar_t yo, scalar_t w)
    {
      SCITBX_ASSERT(grad_yc.size() == n_params);
      SCITBX_ASSERT(!finalised());
      add_equation(yc, grad_yc.begin(), yo, w);
    }

    /// Overload for when efficiency is paramount.
    /** This shall not be called after finalise() has been called
        but this is not enforced for speed.
     */
    void add_equation(scalar_t yc, scalar_t const *grad_yc,
                      scalar_t yo, scalar_t w)
    {
      add_residual(yc, yo, w);
      grad_yc_dot_grad_yc.add(grad_yc, w);
      for (int i=0; i<n_params; ++i) {
        yo_dot_grad_yc[i] += w * yo * grad_yc[i];
        yc_dot_grad_yc[i] += w * yc * grad_yc[i];
      }
    }

    /** \brief Add an equation whose gradients are written here rather than
        handed over.

    open_equation() returns n_parameters() scalars for the caller to fill with
    \f$\nabla y_c\f$, and commit_equation() then does exactly what add_equation()
    does with a vector it was given. Splitting it this way lets the caller build
    the gradients in place: the accumulator has to have them in a row of its own
    whatever happens, so a caller which can write them there directly saves a
    copy of the vector per equation.

    The arithmetic is the one add_equation() does, in the order it does it --
    the vector sums read the unweighted gradients, and the row is weighted after
    -- so the two give the same normal equations to the last bit.
     */
    //@{
    scalar_t *open_equation() {
      return grad_yc_dot_grad_yc.open_row();
    }
    void commit_equation(scalar_t yc, scalar_t const *grad_yc,
                         scalar_t yo, scalar_t w)
    {
      add_residual(yc, yo, w);
      for (int i=0; i<n_params; ++i) {
        yo_dot_grad_yc[i] += w * yo * grad_yc[i];
        yc_dot_grad_yc[i] += w * yc * grad_yc[i];
      }
      grad_yc_dot_grad_yc.commit_row(w);
    }
    //@}

    /// Add many equations in one go
    void add_equations(af::const_ref<scalar_t> const &yc,
                       af::const_ref<scalar_t, af::mat_grid> const &jacobian_yc,
                       af::const_ref<scalar_t> const &yo,
                       af::const_ref<scalar_t> const &w)
    {
      SCITBX_ASSERT(   yc.size() == jacobian_yc.n_rows()
                    && (!w.size() || yc.size() == w.size()))
                   (yc.size())(jacobian_yc.n_rows())(w.size());
      SCITBX_ASSERT(jacobian_yc.n_columns() == n_parameters())
                   (jacobian_yc.n_columns())(n_parameters());
      for (int i=0; i<yc.size(); ++i) {
        add_equation(yc[i], &jacobian_yc(i, 0), yo[i], w.size() ? w[i] : 1);
      }
    }

    /// Per-thread scratch the OpenMP accumulation wants for the normal matrix
    /** Zero for an accumulator which folds the rows in for itself, and the
        caller should not allocate what it will not use: on a large problem
        those per-thread matrices come to more than everything else together.
     */
    static std::size_t omp_matrix_scratch(int n_parameters, int threads) {
      return sum_of_rank_1_updates_t::omp_matrix_scratch(n_parameters, threads);
    }

#if defined(_OPENMP)
    /* These two touch nothing but the scalar sums, so they are the same
    whichever accumulator is in use. add_equations_omp does touch the matrix and
    is specialised per accumulator in normal_equations_omp.h.
    */
    void add_residuals_omp(const int& n,
      const int& start,
      const int& threads,
      af::const_ref<scalar_t> const& yc,
      af::const_ref<scalar_t> const& yo,
      af::const_ref<scalar_t> const& w)
    {
      n_data += n;
      scalar_t temp2 = 0, temp3 = 0, temp4 = 0;
#pragma omp parallel for reduction(+:temp2, temp3, temp4) num_threads(threads)
      for (int i = start; i < start + n; i++) {
        scalar_t const temp1 = w[i] * yo[i];
        temp2 += temp1 * yo[i];
        temp3 += temp1 * yc[i];
        temp4 += w[i] * yc[i] * yc[i];
      }
      yo_sq += temp2;
      yo_dot_yc += temp3;
      yc_sq += temp4;
    }

    void add_residuals_omp(const int& n,
      const int& start,
      const int& threads,
      af::const_ref<scalar_t> const& yc,
      af::const_ref<scalar_t> const& yo)
    {
      n_data += n;
      scalar_t temp1 = 0, temp2 = 0, temp3 = 0;
#pragma omp parallel for reduction(+:temp1, temp2, temp3) num_threads(threads)
      for (int i = start; i < start + n; i++) {
        temp1 += yo[i] * yo[i];
        temp2 += yo[i] * yc[i];
        temp3 += yc[i] * yc[i];
      }
      yo_sq += temp1;
      yo_dot_yc += temp2;
      yc_sq += temp3;
    }

    /// Add many equations in one go using OpenMP
    void add_equations_omp(const int& n_ref,
      const int& n_par,
      const int& chunk_size,
      const int& start,
      const int& threads,
      std::vector<FloatType>& matrix,          //
      std::vector<FloatType>& yo_dot_grad_yc_, // These are arrays passed as locals for each thread to reduce overhead of generating them for each chunk
      std::vector<FloatType>& yc_dot_grad_yc_, //
      af::const_ref<scalar_t> const& yc,
      std::vector<FloatType> const& jacobian_yc,
      af::const_ref<scalar_t> const& yo,
      af::const_ref<scalar_t> const& w)
    {
      throw SCITBX_NOT_IMPLEMENTED();
    }
#endif

    /// Addition in the sense of the L.S. objective functions
    /**
     *  The overall factor for this objective and other's objective are the same.
     *  In term of normal equations, this appends other's equations to this.
     */

    non_linear_ls_with_separable_scale_factor
    &operator+=(non_linear_ls_with_separable_scale_factor const &other) {
      SCITBX_ASSERT(!finalised());
      SCITBX_ASSERT(!other.finalised());
      n_data += other.n_data;
      yo_sq += other.yo_sq;
      yo_dot_yc += other.yo_dot_yc;
      yc_sq += other.yc_sq;

      grad_yc_dot_grad_yc += other.grad_yc_dot_grad_yc;
      yo_dot_grad_yc += other.yo_dot_grad_yc;
      yc_dot_grad_yc += other.yc_dot_grad_yc;

      return *this;
    }

    /// \f$\sum w y_o^2\f$
    /** This is the normalisation that guarantees
        that \f$L(K, x)\f$ is between 0 and 1.
     */
    scalar_t sum_w_yo_sq() const {
      SCITBX_ASSERT(finalised());
      return yo_sq;
    }

    /** \brief The value \f$ K^*(x) \f$ of the scale factor optimising the L.S. objective for a given constant \f$ x \f$.
     */
    scalar_t optimal_scale_factor() const {
      SCITBX_ASSERT(finalised());
      return yo_dot_yc/yc_sq;
    }

    /// The value of the minimised function, for the optimal scale factor
    /** This is \f$L(K^*(x), x)\f$, plus the contributions added to
        the reduced_problem().
     */
    scalar_t objective() const {
      SCITBX_ASSERT(finalised());
      return reduced_ls.objective();
    }

    /// \f$\chi^2\f$ of the fit.
    /** The \f$\chi^2\f$ for the fit of \f$y_c(x)\f$ to \f$y_o\f$.
        This does include the contributions added to the reduced_problem().
     */
    scalar_t chi_sq() const {
      SCITBX_ASSERT(finalised());
      return (r_sq + 2*(reduced_ls.objective() - objective_))/dof();
    }

    /// Equation accumulation is finished.
    /** The reduced normal equations for \f$ x \f$ as per step 2 are constructed
     */
    void finalise(bool objective_only=false) {
      SCITBX_ASSERT(!finalised() && n_equations())(n_equations());
      finalised_ = true;

      grad_yc_dot_grad_yc.finalise();
      a = grad_yc_dot_grad_yc;

      scalar_t k_star = optimal_scale_factor(), k_star_sq = k_star*k_star;
      r_sq = yo_sq*(1 - (k_star_sq * yc_sq)/yo_sq);
      objective_ = r_sq/2;
      if (normalised()) objective_ /= yo_sq;

      vector_owning_ref_t b = yo_dot_grad_yc;
      reduced_ls = non_linear_ls<scalar_t>(n_data,
                                           objective_, b.array(), a.array());

      if (objective_only) return;

      scalar_t r_dot_yc = yo_dot_yc - k_star*yc_sq;
      scalar_t inv_yc_sq = 1./yc_sq;
      for (int i=0; i<n_params; ++i) {
        scalar_t r_dot_grad_yc_i = yo_dot_grad_yc[i] - k_star*yc_dot_grad_yc[i];
        grad_k_star[i] = inv_yc_sq*(r_dot_grad_yc_i - k_star*yc_dot_grad_yc[i]);
        b[i] = k_star*r_dot_grad_yc_i + grad_k_star[i]*r_dot_yc;
      }
      scalar_t *pa = a.begin();
      for (int i=0; i<n_params; ++i) for (int j=i; j<n_params; ++j) {
        scalar_t a_ij = *pa;
        a_ij = k_star_sq*a_ij
             + k_star*(  yc_dot_grad_yc[i]*grad_k_star[j]
                       + yc_dot_grad_yc[j]*grad_k_star[i])
             + grad_k_star[i]*grad_k_star[j]*yc_sq;
        *pa++ = a_ij;
      }
      if (normalised()) {
        a /= yo_sq;
        b /= yo_sq;
      }
    }

    /// Whether finalise has been called.
    bool finalised() const { return finalised_; }

    /// The linear L.S. problem to solve for a step toward the minimum.
    linear_ls<scalar_t> &step_equations() {
      return reduced_problem().step_equations();
    }

    /// The non-linear problem with the scale factor already optimised away
    /** The main use of this function comes for an objective function

        \[ \tilde{L}(K, x) = L(K, x) + \frac{1}{2} \|r(x\|^2 \]

        for some residual vector \f$r(x)\f$ independent of the overall scale
        factor that the first term depends upon. The equations for that
        second term may then be accumulated into the object returned by this
        function, to produce the correct equations for
        \f$\tilde{L}(K^*(x), x)\f$.

        This would not be possible with step_equations() which looses
        sight of the objective value.

        Invariant: reduced_problem().step_equations() and this->step_equations()
        are identical (i.e. modify one, modifies the other).
     */
    non_linear_ls<scalar_t> &reduced_problem() {
      SCITBX_ASSERT(finalised());
      return reduced_ls;
    }

    /// Ready this for another computation of the normal equations
    void reset() {
      n_data = 0;
      yo_dot_yc = 0; yc_sq = 0; yo_sq = 0;
      grad_yc_dot_grad_yc.reset();
      std::fill(yo_dot_grad_yc.begin(), yo_dot_grad_yc.end(), scalar_t(0));
      std::fill(yc_dot_grad_yc.begin(), yc_dot_grad_yc.end(), scalar_t(0));
      std::fill(grad_k_star.begin(), grad_k_star.end(), scalar_t(0));
      finalised_ = false;
    }

  private:
    scalar_t yo_dot_yc, yo_sq, yc_sq, r_sq, objective_;
    int n_params;
    std::size_t n_data;
    bool normalised_;
    sum_of_rank_1_updates_t grad_yc_dot_grad_yc;
    symmetric_matrix_owning_ref_t a; // normal matrix stored
                                     // as packed upper diagonal
    vector_owning_ref_t yo_dot_grad_yc, yc_dot_grad_yc, grad_k_star;
    bool finalised_;
    non_linear_ls<scalar_t> reduced_ls;
  };

  /// Non-linear L.S. with a scale factor that is given, not optimised away.
  /** Same model \f$y_o \approx K y_c(x)\f$ and the same accumulation as
      non_linear_ls_with_separable_scale_factor, but \f$K\f$ is held fixed:

      \f[ A = K^2 \sum_i w_i \nabla y_{c,i} \nabla y_{c,i}^T, \quad
          b = K \sum_i w_i (y_{o,i} - K y_{c,i}) \nabla y_{c,i} \f]

      This exists for maximum-likelihood refinement. There the accumulator is
      fed an *effective* observation and weight derived from the likelihood, so
      that \f$ w(y_o - K y_c) \f$ is the likelihood gradient; variable
      projection would then solve for a second scale factor on top of the one
      already inside the likelihood's alpha and beta, and the fixed point of
      that is not the maximum of the likelihood.

      The public surface is deliberately identical to the separable class, down
      to optimal_scale_factor(), so that everything built on that one -
      smtbx's crystallographic_ls, the restraints, the solvers - takes this by
      substitution rather than by special-casing.
   */
  template <typename FloatType,
            template<typename> class SumOfRank1Updates=matrix::sum_of_symmetric_rank_1_updates>
  class non_linear_ls_with_fixed_scale_factor
  {
  public:
    SCITBX_LSTBX_DECLARE_ARRAY_TYPE(FloatType);
    typedef SumOfRank1Updates<FloatType> sum_of_rank_1_updates_t;

    non_linear_ls_with_fixed_scale_factor(
      int n_parameters,
      bool normalised=false,
      std::size_t accumulator_buffer_bytes=0)
      : yo_dot_yc(0), yc_sq(0), yo_sq(0),
        k_(1),
        n_params(n_parameters),
        n_data(0),
        normalised_(normalised),
        grad_yc_dot_grad_yc(n_parameters, accumulator_buffer_bytes),
        yo_dot_grad_yc(n_parameters),
        yc_dot_grad_yc(n_parameters),
        finalised_(false),
        reduced_ls(n_parameters)
    {}

    int n_parameters() const { return n_params; }

    std::size_t n_equations() const {
      return finalised() ? reduced_ls.n_equations() : n_data;
    }

    std::size_t dof() const { return n_equations() - n_parameters(); }

    bool normalised() const { return normalised_; }

    /// The scale factor this works at. Set it before finalise().
    void set_scale_factor(scalar_t k) {
      SCITBX_ASSERT(!finalised());
      k_ = k;
    }

    void add_residual(scalar_t yc, scalar_t yo, scalar_t w) {
      n_data++;
      yo_sq += w * yo * yo;
      yo_dot_yc += w * yo * yc;
      yc_sq += w * yc * yc;
    }

    void add_equation(scalar_t yc, af::const_ref<scalar_t> const &grad_yc,
                      scalar_t yo, scalar_t w)
    {
      SCITBX_ASSERT(grad_yc.size() == n_params);
      SCITBX_ASSERT(!finalised());
      add_equation(yc, grad_yc.begin(), yo, w);
    }

    void add_equation(scalar_t yc, scalar_t const *grad_yc,
                      scalar_t yo, scalar_t w)
    {
      add_residual(yc, yo, w);
      grad_yc_dot_grad_yc.add(grad_yc, w);
      for (int i=0; i<n_params; ++i) {
        yo_dot_grad_yc[i] += w * yo * grad_yc[i];
        yc_dot_grad_yc[i] += w * yc * grad_yc[i];
      }
    }

    scalar_t *open_equation() {
      return grad_yc_dot_grad_yc.open_row();
    }
    void commit_equation(scalar_t yc, scalar_t const *grad_yc,
                         scalar_t yo, scalar_t w)
    {
      add_residual(yc, yo, w);
      for (int i=0; i<n_params; ++i) {
        yo_dot_grad_yc[i] += w * yo * grad_yc[i];
        yc_dot_grad_yc[i] += w * yc * grad_yc[i];
      }
      grad_yc_dot_grad_yc.commit_row(w);
    }

    void add_equations(af::const_ref<scalar_t> const &yc,
                       af::const_ref<scalar_t, af::mat_grid> const &jacobian_yc,
                       af::const_ref<scalar_t> const &yo,
                       af::const_ref<scalar_t> const &w)
    {
      SCITBX_ASSERT(   yc.size() == jacobian_yc.n_rows()
                    && (!w.size() || yc.size() == w.size()))
                   (yc.size())(jacobian_yc.n_rows())(w.size());
      SCITBX_ASSERT(jacobian_yc.n_columns() == n_parameters())
                   (jacobian_yc.n_columns())(n_parameters());
      for (int i=0; i<yc.size(); ++i) {
        add_equation(yc[i], &jacobian_yc(i, 0), yo[i], w.size() ? w[i] : 1);
      }
    }

    static std::size_t omp_matrix_scratch(int n_parameters, int threads) {
      return sum_of_rank_1_updates_t::omp_matrix_scratch(n_parameters, threads);
    }

#if defined(_OPENMP)
    void add_residuals_omp(const int& n,
      const int& start,
      const int& threads,
      af::const_ref<scalar_t> const& yc,
      af::const_ref<scalar_t> const& yo,
      af::const_ref<scalar_t> const& w)
    {
      n_data += n;
      scalar_t temp2 = 0, temp3 = 0, temp4 = 0;
#pragma omp parallel for reduction(+:temp2, temp3, temp4) num_threads(threads)
      for (int i = start; i < start + n; i++) {
        scalar_t const temp1 = w[i] * yo[i];
        temp2 += temp1 * yo[i];
        temp3 += temp1 * yc[i];
        temp4 += w[i] * yc[i] * yc[i];
      }
      yo_sq += temp2;
      yo_dot_yc += temp3;
      yc_sq += temp4;
    }

    void add_residuals_omp(const int& n,
      const int& start,
      const int& threads,
      af::const_ref<scalar_t> const& yc,
      af::const_ref<scalar_t> const& yo)
    {
      n_data += n;
      scalar_t temp1 = 0, temp2 = 0, temp3 = 0;
#pragma omp parallel for reduction(+:temp1, temp2, temp3) num_threads(threads)
      for (int i = start; i < start + n; i++) {
        temp1 += yo[i] * yo[i];
        temp2 += yo[i] * yc[i];
        temp3 += yc[i] * yc[i];
      }
      yo_sq += temp1;
      yo_dot_yc += temp2;
      yc_sq += temp3;
    }

    /// Not provided: the caller must not take the OpenMP path with this.
    /** The separable class has this specialised per accumulator in
        normal_equations_omp.h. Rather than a second set of specialisations for
        a path maximum likelihood does not need yet, this refuses loudly - which
        is better than a silently wrong normal matrix.
     */
    void add_equations_omp(const int& n_ref,
      const int& n_par,
      const int& chunk_size,
      const int& start,
      const int& threads,
      std::vector<FloatType>& matrix,
      std::vector<FloatType>& yo_dot_grad_yc_,
      std::vector<FloatType>& yc_dot_grad_yc_,
      af::const_ref<scalar_t> const& yc,
      std::vector<FloatType> const& jacobian_yc,
      af::const_ref<scalar_t> const& yo,
      af::const_ref<scalar_t> const& w)
    {
      throw SCITBX_NOT_IMPLEMENTED();
    }
#endif

    non_linear_ls_with_fixed_scale_factor
    &operator+=(non_linear_ls_with_fixed_scale_factor const &other) {
      SCITBX_ASSERT(!finalised());
      SCITBX_ASSERT(!other.finalised());
      n_data += other.n_data;
      yo_sq += other.yo_sq;
      yo_dot_yc += other.yo_dot_yc;
      yc_sq += other.yc_sq;

      grad_yc_dot_grad_yc += other.grad_yc_dot_grad_yc;
      yo_dot_grad_yc += other.yo_dot_grad_yc;
      yc_dot_grad_yc += other.yc_dot_grad_yc;

      return *this;
    }

    scalar_t sum_w_yo_sq() const {
      SCITBX_ASSERT(finalised());
      return yo_sq;
    }

    /// The scale factor in use. Fixed, so nothing was optimised to get it.
    /** Named as in the separable class on purpose: callers ask for the scale
        the equations were built at, and that question has an answer here too.
     */
    scalar_t optimal_scale_factor() const {
      return k_;
    }

    scalar_t objective() const {
      SCITBX_ASSERT(finalised());
      return reduced_ls.objective();
    }

    scalar_t chi_sq() const {
      SCITBX_ASSERT(finalised());
      return (r_sq + 2*(reduced_ls.objective() - objective_))/dof();
    }

    void finalise(bool objective_only=false) {
      SCITBX_ASSERT(!finalised() && n_equations())(n_equations());
      finalised_ = true;

      grad_yc_dot_grad_yc.finalise();
      a = grad_yc_dot_grad_yc;

      scalar_t k_sq = k_*k_;
      r_sq = yo_sq - 2*k_*yo_dot_yc + k_sq*yc_sq;
      objective_ = r_sq/2;
      if (normalised()) objective_ /= yo_sq;

      vector_owning_ref_t b = yo_dot_grad_yc;
      reduced_ls = non_linear_ls<scalar_t>(n_data,
                                           objective_, b.array(), a.array());

      if (objective_only) return;

      for (int i=0; i<n_params; ++i) {
        b[i] = k_*(yo_dot_grad_yc[i] - k_*yc_dot_grad_yc[i]);
      }
      if (k_ != 1) {
        scalar_t *pa = a.begin();
        for (int i=0; i<n_params; ++i) for (int j=i; j<n_params; ++j) {
          *pa++ *= k_sq;
        }
      }
      if (normalised()) {
        a /= yo_sq;
        b /= yo_sq;
      }
    }

    bool finalised() const { return finalised_; }

    linear_ls<scalar_t> &step_equations() {
      return reduced_problem().step_equations();
    }

    non_linear_ls<scalar_t> &reduced_problem() {
      SCITBX_ASSERT(finalised());
      return reduced_ls;
    }

    void reset() {
      n_data = 0;
      yo_dot_yc = 0; yc_sq = 0; yo_sq = 0;
      grad_yc_dot_grad_yc.reset();
      std::fill(yo_dot_grad_yc.begin(), yo_dot_grad_yc.end(), scalar_t(0));
      std::fill(yc_dot_grad_yc.begin(), yc_dot_grad_yc.end(), scalar_t(0));
      finalised_ = false;
    }

  private:
    scalar_t yo_dot_yc, yo_sq, yc_sq, r_sq, objective_, k_;
    int n_params;
    std::size_t n_data;
    bool normalised_;
    sum_of_rank_1_updates_t grad_yc_dot_grad_yc;
    symmetric_matrix_owning_ref_t a;
    vector_owning_ref_t yo_dot_grad_yc, yc_dot_grad_yc;
    bool finalised_;
    non_linear_ls<scalar_t> reduced_ls;
  };

// OpenMP specialisation
#if defined(_OPENMP)
  #include <scitbx/lstbx/normal_equations_omp.h>
#endif
}}}

#endif // GUARD
