/// Tools implementing the Gauss-Newton method for non-linear least-squares.

#ifndef SCITBX_GAUSS_NEWTON_H
#define SCITBX_GAUSS_NEWTON_H

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/shared_algebra.h>
#include <scitbx/array_family/ref_algebra.h>
#include <scitbx/array_family/owning_ref.h>
#include <scitbx/array_family/accessors/row_and_column.h>
#include <scitbx/matrix/cholesky.h>
#include <scitbx/matrix/matrix_vector_operations.h>
#include <scitbx/sparse/matrix.h>
#include <scitbx/sparse/triangular.h>

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
        normal_matrix_(n_parameters),
        right_hand_side_(n_parameters)
    {}

    /// Number of unknown parameters
    int n_parameters() { return right_hand_side_.size(); }

    /// Initialise the least-squares problem with the given normal matrix A
    /// and right hand side b
    linear_ls(symmetric_matrix_t const &a, vector_t const &b)
      : solved_(false),
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
      scalar_t *p = normal_matrix_.begin();
      for (int i=0; i<n_parameters(); ++i)  {
        right_hand_side_[i] += w * a_row[i] * b_i;
        for (int j=i; j<n_parameters(); ++j) *p++ += w * a_row[i] * a_row[j];
      }
    }

    /// Add the equations A x = b with the given weights
    /** w[i] weights the i-th equation, i.e. the row \f$ A_{i.} \f$.
        If negate_right_hand_side, then the equation is A x + b = 0 instead
     */
    void add_equations(af::const_ref<scalar_t> const &b,
                       sparse::matrix<scalar_t> const &a,
                       af::const_ref<scalar_t> const &w,
                       bool negate_right_hand_side=false)
    {
      SCITBX_ASSERT(   a.n_rows() == b.size()
                    && b.size()   == w.size())(a.n_rows())(b.size())(w.size());
      SCITBX_ASSERT(a.n_cols() == n_parameters());
      sparse::matrix<scalar_t>
      at_w_a = a.this_transpose_times_diagonal_times_this(w);
      vector_t a_t_w_b = a.transpose_times((w * b).const_ref());
      normal_matrix_ += sparse::upper_diagonal_of(at_w_a);
      if (negate_right_hand_side) right_hand_side_ -= a_t_w_b.const_ref();
      else                        right_hand_side_ += a_t_w_b.const_ref();

    }

    /// Reset the state to construction time, i.e. no equations accumulated
    void reset() {
      solved_ = false;
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
      u_transpose_u_decomposition_in_place<scalar_t> chol(normal_matrix_);
      SCITBX_ASSERT(!chol.failure);
      chol.solve_in_place(right_hand_side_);
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

  private:
    bool solved_;
    symmetric_matrix_owning_ref_t normal_matrix_;
    vector_owning_ref_t right_hand_side_;
  };


  /// Normal equations for non-linear least-squares
  /** The least-squares target reads

      \f[ L(x) = \|r_i(x)\|^2 \f]

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
      : linearised(n_parameters),
        r_sq(0)
    {}

    /// Construct with an exiting L.S. problem.
    /** That is
          - this->objective()  == objective
          - this->step_equations().right_hand_side() == opposite_of_grad_objective
          - this->step_equations().normal_matrix() == normal_matrix
     */
    non_linear_ls(scalar_t objective,
                  vector_t const &opposite_of_grad_objective,
                  symmetric_matrix_t const &normal_matrix)
    : r_sq(2*objective), linearised(normal_matrix, opposite_of_grad_objective)
    {}

    /// Number of unknown parameters
    int n_parameters() { return linearised.n_parameters(); }

    /// Add the given residual with the given weight
    void add_residual(scalar_t r, scalar_t w) {
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
                       af::const_ref<scalar_t> w)
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

    /// Objective value \f$L(x)\f$ for the current value of the unknowns
    scalar_t objective() { return r_sq/2; }

    /// Linearised equations to solve for a step
    linear_ls<scalar_t> &step_equations() { return linearised; }

    /// Reset the state to construction time, i.e. no equations accumulated
    void reset() {
      linearised.reset();
      r_sq = 0;
    }

  private:
    scalar_t r_sq;
    linear_ls<scalar_t> linearised;
  };


  /// Normal equations for least-squares fit with an overall scale.
  /** The least-squares target reads

      \f[ L(K, x) = \frac{ \sum w ( K y_c(x) - y_o )^2 }{ \sum w y_o^2 } \f]

      where the both of \f$ y_c(x) \f$ and \f$ y_o \f$ are vectors,
      respectively the model to fit to the data. Alternatively, the
      non-normalised

      \f[ \tilde{L}(K, x) = \sum w ( K y_c(x) - y_o )^2 \f]

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
  template <typename FloatType>
  class non_linear_ls_with_separable_scale_factor
  {
  public:
    SCITBX_LSTBX_DECLARE_ARRAY_TYPE(FloatType);

    /// Construct a least-squares problem with the given number of parameters.
    /** That is the length of the vector \f$ x \f$. The flag normalised
        specify whether to use the normalised objective \f$L\f$ or the
        non-normalised objective \f$\tilde{L}\f$.
     */
    non_linear_ls_with_separable_scale_factor(int n_parameters,
                                              bool normalised=true)
      : yo_dot_yc(0), yc_sq(0), yo_sq(0),
        n_params(n_parameters),
        normalised_(normalised),
        a(n_parameters),
        yo_dot_grad_yc(n_parameters),
        yc_dot_grad_yc(n_parameters),
        grad_k_star(n_parameters),
        finalised_(false),
        reduced_ls(n_parameters)
    {}

    /// Number of unknown parameters, not including the overall scale factor
    int n_parameters() { return n_params; }

    /// Whether the L.S. target is normalised by \f$ \sum w y_o^2 \f$ or not
    bool normalised() { return normalised_; }

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
      yo_sq += w * yo * yo;
      yo_dot_yc += w * yo * yc;
      yc_sq += w * yc * yc;
      double *pa = a.begin();
      for (int i=0; i<n_params; ++i) for (int j=i; j<n_params; ++j) {
        *pa++ += w * grad_yc[i] * grad_yc[j];
      }
      for (int i=0; i<n_params; ++i) {
        yo_dot_grad_yc[i] += w * yo * grad_yc[i];
        yc_dot_grad_yc[i] += w * yc * grad_yc[i];
      }
    }

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

    /// \sum w y_o^2
    /** This is the normalisation that guarantees
        that \f$L(K, x)\f$ is between 0 and 1.
     */
    scalar_t sum_w_yo_sq() {
      SCITBX_ASSERT(finalised());
      return yo_sq;
    }

    /** \brief The value \f$ K^*(x) \f$ of the scale factor optimising the L.S. objective for a given constant \f$ x \f$.
     */
    scalar_t optimal_scale_factor() {
      SCITBX_ASSERT(finalised());
      return yo_dot_yc/yc_sq;
    }

    /** \brief The value of \f$L(K, x)\f$
     for the optimised scale factor \f$ K^*(x) \f$ and the input \f$yc_(x)\f$.
     */
    scalar_t objective() {
      return reduced_problem().objective();
    }

    /// Equation accumulation is finished.
    /** The reduced normal equations for \f$ x \f$ as per step 2 are constructed
     */
    void finalise() {
      SCITBX_ASSERT(!finalised());
      finalised_ = true;

      scalar_t k_star = optimal_scale_factor(), k_star_sq = k_star*k_star;
      scalar_t objective_ = (1 - (k_star_sq * yc_sq)/yo_sq)/2;
      if (!normalised()) objective_ *= yo_sq;

      vector_owning_ref_t b = yo_dot_grad_yc;
      reduced_ls = non_linear_ls<scalar_t>(objective_, b.array(), a.array());

      scalar_t r_dot_yc = yo_dot_yc - k_star*yc_sq;
      scalar_t inv_yc_sq = 1./yc_sq;
      for (int i=0; i<n_params; ++i) {
        scalar_t r_dot_grad_yc_i = yo_dot_grad_yc[i] - k_star*yc_dot_grad_yc[i];
        grad_k_star[i] = inv_yc_sq*(r_dot_grad_yc_i - k_star*yc_dot_grad_yc[i]);
        b[i] = k_star*r_dot_grad_yc_i + grad_k_star[i]*r_dot_yc;
      }
      double *pa = a.begin();
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
    bool finalised() {
      return finalised_;
    }

    /// Reduced normal equations
    linear_ls<scalar_t> &step_equations() {
      return reduced_problem().step_equations();
    }

    /// The non-linear problem with the scale factor already optimised away
    /** The main use of this function comes for an objective function

        \[ \tilde{L}(K, x) = L(K, x) + \|r(x\|^2 \]

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
      yo_dot_yc = 0; yc_sq = 0; yo_sq = 0;
      std::fill(a.begin(), a.end(), scalar_t(0));
      std::fill(yo_dot_grad_yc.begin(), yo_dot_grad_yc.end(), scalar_t(0));
      std::fill(yc_dot_grad_yc.begin(), yc_dot_grad_yc.end(), scalar_t(0));
      std::fill(grad_k_star.begin(), grad_k_star.end(), scalar_t(0));
      finalised_ = false;
    }

  private:
    scalar_t yo_dot_yc, yo_sq, yc_sq;
    int n_params;
    bool normalised_;
    symmetric_matrix_owning_ref_t a; // normal matrix stored
                                     // as packed upper diagonal
    vector_owning_ref_t yo_dot_grad_yc, yc_dot_grad_yc, grad_k_star;
    bool finalised_;
    non_linear_ls<scalar_t> reduced_ls;
  };



  /// Normal equations for least-squares with a separable linear part
  /** The least-squares target reads

   \f$ L(a, x) = \sum w ( y_o - Y_c(x) a )

   Thus the model \f$ y_c \f$ fitted to the data \f$ y_o \f$ is the linear
   combination of some \f$ y_{c,1}, \cdots, y_{c,p} \f$ with coefficients
   \f$ a_1, \cdots, a_p \f$.

   The separability of the problem is exploited by

   - first computing \f$ a^* = \argmin_a L(a, x) \f$,

   - and then building the Newton equations for \f$ L(a^*(x), x) \f$
     in the Gauss approximation of small residuals.

   TODO: fix implementation as normal_equations_separating_scale_factor
   has been fixed.
   */
  template <typename FloatType>
  class separable_non_linear_ls
  {
  public:
    SCITBX_LSTBX_DECLARE_ARRAY_TYPE(FloatType);

    separable_non_linear_ls(int n_parameters, int n_bases)
      : n(n_parameters), p(n_bases),
        linear_part(p),
        grad_yc_trans_grad_yc(n*(n+1)/2 * p*p),
        yo_trans_grad_yc(af::mat_grid(n,p)),
        yc_trans_grad_yc(n*p*p),
        a(n), b(n),
        has_solved_separable_linear_part(false)
    {}

    void add_equation(af::const_ref<scalar_t> const &yc,
                      af::const_ref<scalar_t, af::mat_grid> const &grad_yc,
                      scalar_t yo, scalar_t w)
    {
      linear_part.add_equation(yo, yc, w);

      double *q = grad_yc_trans_grad_yc.begin();
      for (int i=0; i<n; ++i) for (int j=i; j<n; ++j) {
        for (int k=0; k<p; ++k) for (int l=0; l<p; ++l) {
          *q++ += w * grad_yc(i,k) * grad_yc(j,l);
        }
      }

      double *yo_tr_grad_yc = yo_trans_grad_yc.begin();
      double *yc_tr_grad_yc = yc_trans_grad_yc.begin();
      for (int i=0; i<n; ++i) {
        for (int k=0; k<p; ++k) {
          *yo_tr_grad_yc++ += w * yo * grad_yc(i,k);
          for (int l=0; l<p; ++l) {
            *yc_tr_grad_yc++ += w * yc(k) * grad_yc(i,l);
          }
        }
      }
    }

    vector_t optimal_linear_coefficients() {
      if(!has_solved_separable_linear_part) optimise_separable_linear_part();
      return a_star.array();
    }

    linear_ls<scalar_t> equations() {
      using namespace scitbx::matrix;
      if(!has_solved_separable_linear_part) optimise_separable_linear_part();

      int sym_p = p*(p+1)/2;
      scalar_t *a_star_ = a_star.begin();
      scalar_t *alpha0 = yo_trans_grad_yc.begin();

      scalar_t *alpha = alpha0;
      scalar_t *yc_tr_grad_yc = yc_trans_grad_yc.begin();
      for (int i=0; i<n; ++i, yc_tr_grad_yc += p*p, alpha += p)
      {
        scalar_t *beta = alpha;
        matrix_vector(p, p, yc_tr_grad_yc, a_star_, beta, -1., 1.);
        b[i] = dot(p, a_star_, beta);
        matrix_transposed_vector(p, p, yc_tr_grad_yc, a_star_, alpha, -1., 1.);
        forward_substitution_given_transpose(p, lin_u.begin(), alpha);
      }

      scalar_t *a_ = a.begin();
      alpha = alpha0;
      scalar_t *grad_yc_tr_grad_yc = grad_yc_trans_grad_yc.begin();
      for (int i=0; i<n; ++i) for (int j=i; j<n; ++j) {
        scalar_t a_ij = quadratic_form(p, grad_yc_tr_grad_yc, a_star_)
                        -dot(p, &alpha[i*p], &alpha[j*p]);
        *a_++ = a_ij;
        grad_yc_tr_grad_yc += p*p;
      }

      return linear_ls<scalar_t>(a.array(), b.array());
    }

  private:
    void optimise_separable_linear_part() {
      using namespace scitbx::matrix;
      symmetric_matrix_owning_ref_t lin_a = linear_part.normal_matrix();
      vector_owning_ref_t lin_b = linear_part.right_hand_side();
      cholesky::u_transpose_u_decomposition_in_place<scalar_t> chol(lin_a);
      SCITBX_ASSERT(!chol.failure);
      lin_u = lin_a.array();
      cholesky::solve_in_place::using_u_transpose_u(lin_u, lin_b);
      a_star = lin_b.array();
      has_solved_separable_linear_part = true;
    }

    int n;
    int p;

    linear_ls<scalar_t> linear_part;

    af::ref_owning_shared<scalar_t> grad_yc_trans_grad_yc;
    matrix_owning_ref_t yo_trans_grad_yc;
    af::ref_owning_shared<scalar_t> yc_trans_grad_yc;

    vector_owning_ref_t b;
    symmetric_matrix_owning_ref_t a;

    bool has_solved_separable_linear_part;
    symmetric_matrix_owning_ref_t lin_u;
    vector_owning_ref_t a_star;
  };



}}}

#endif // GUARD
