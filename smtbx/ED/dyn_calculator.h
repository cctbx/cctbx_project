#pragma once
#include <smtbx/ED/ed.h>
#include <smtbx/ED/hermitian_eigen.h>
#include <boost/shared_ptr.hpp>
namespace smtbx { namespace ED
{
  using namespace cctbx;

  /** @brief Dynamical diffraction: what the calculation is, before the code.

  An electron passing through a crystal is scattered many times over, so the
  intensity of a reflection is not proportional to |F(h)|^2 as it is for X-rays.
  The beams exchange amplitude with each other all the way through the crystal,
  and what leaves the far side is the solution of a coupled system rather than
  a single Fourier coefficient.

  Take the n beams thought to matter for one reflection -- the incident beam
  first, then the others -- and let psi(t) hold their complex amplitudes at
  depth t. To a good approximation they obey

      d psi / dt = i (pi / K) A psi

  with A the scattering matrix: off the diagonal A(i,j) is the Fourier
  coefficient of the crystal potential connecting beams i and j, which is the
  kinematic structure factor of h_i - h_j, and on the diagonal it carries how
  far beam i sits from its exact Bragg condition (its excitation error). A is
  Hermitian, which is the whole reason the arithmetic below is affordable.
  See utils::build_Ug_matrix and each subclass's build() for the two halves.

  A does not depend on t, so the system integrates in closed form:

      psi(t) = exp(i pi t A / K) psi(0),   psi(0) = (1, 0, ... 0)

  -- everything starts in the incident beam. So the amplitude of beam j is
  element (j,0) of a matrix exponential, and the measured intensity is its
  modulus squared. Diagonalising the Hermitian A once, A = U L U*, turns the
  exponential into a scalar one per eigenvalue:

      exp(i pi t A / K) = U diag(exp(i pi t L_k / K)) U*

  which is what every calc_amps* here does: one call to hermitian_eigen, then a
  weighted sum over eigenvectors. (It overwrites A with U, so `A` means the
  scattering matrix before that call and the eigenvectors after it.)

  Refinement needs the derivative of the intensity with respect to each refined
  parameter, and that is where it stops being a textbook exercise. Moving an
  atom changes every structure factor, hence every element of A, so what is
  needed is the derivative of a matrix exponential with respect to a matrix
  that does not commute with itself at different parameters -- for which
  d/dp exp(cA) = c exp(cA) A' is simply false. The right expression, the
  Frechet derivative in the eigenbasis (Daleckii-Krein; the code cites
  Bernoulli 9(5), 2003, 895-919), is

      d/dp exp(cA) = U [ G o (U* A' U) ] U*

  where A' = dA/dp, `o` is the elementwise product, and G is the matrix of
  divided differences of the scalar exponentials,

      G(i,j) = (exp(c L_i) - exp(c L_j)) / (L_i - L_j),   G(i,i) = c exp(c L_i)

  the diagonal being the limit of the off-diagonal as the eigenvalues meet.
  Every calc_amps_ext* below builds that G. A' for each parameter is assembled
  by utils::build_D_matrices from the kinematic design matrix: element (i,j) of
  A' is d/dp of the structure factor of h_i - h_j, which the kinematic
  calculation has already worked out. Finally, intensity from amplitude:

      dI/dp = d|psi|^2/dp = 2 Re( conj(psi) dpsi/dp )

  which is the last line of each of those functions.
  */
  template <typename FloatType>
  class a_dyn_calculator {
  public:
    ED_UTIL_TYPEDEFS;

    // mat_Ug will be affected if build is called!
    a_dyn_calculator(const af::shared<miller::index<> >& indices,
      const cmat_t& mat_Ug,
      const cart_t& K,
      const mat3_t& RMf,
      const cart_t& N,
      FloatType thickness)
      : indices(indices),
      A(mat_Ug), K(K), RMf(RMf), N(N),
      thickness(thickness)
    {}

    a_dyn_calculator(const af::shared<miller::index<> >& indices,
      const cart_t& K,
      const cart_t& N,
      FloatType thickness)
      : indices(indices),
      K(K), N(N),
      thickness(thickness)
    {}

    virtual ~a_dyn_calculator() {}

    // mat_Ug will be NOT be affected - deep copied
    a_dyn_calculator& reset(const cmat_t& m, const mat3_t& RMf) {
      A = m.deep_copy();
      this->RMf = RMf;
      return build();
    }

    // mat_Ug will be NOT be affected - deep copied
    a_dyn_calculator& reset(const af::shared<miller::index<> > &indices_,
      const cmat_t& m, const mat3_t& RMf)
    {
      A = m.deep_copy();
      indices = indices_;
      this->RMf = RMf;
      return build();
    }

    virtual af::shared<complex_t> calc_amps(size_t num,
      bool include_incident=false) = 0;
    // 0 is for the incident beam
    virtual complex_t calc_amps_1(size_t idx) = 0;

    virtual af::shared<complex_t> calc_amps_ext(
      af::shared<cmat_t> const& Ds_kin,
      bool grad_thickness,
      mat_t& D_dyn,
      size_t num) = 0;

    //D_dyn has one row as output, 0 is the incident beam
    virtual complex_t calc_amps_ext_1(
      const af::shared<cmat_t>& Ds_kin,
      bool grad_thickness,
      mat_t& D_dyn,
      size_t idx) = 0;

    // recomputes the Eigen matrix
    virtual a_dyn_calculator& build() = 0;
    const cart_t K, N;
    const cmat_t& get_matrix() const {
      return A;
    }
  protected:
    /** @brief One output beam's row of dI/dp, for every refined parameter.

    The derivative of the intensity of beam @p idx with respect to parameter p
    is, writing V_p for A* D_p A and G for the divided differences of the
    eigenvalue exponentials,

      dp = sum_ij  A(idx,i) G(i,j) A*(j,0) V_p(i,j)  =  sum_ij W(i,j) V_p(i,j)

    Every factor of W is a property of the eigendecomposition, not of the
    parameter. Substituting V_p and gathering the sums the other way round,

      dp = sum_ab D_p(a,b) sum_ij conj(A(a,i)) W(i,j) A(b,j)
         = sum_ab D_p(a,b) M(a,b),   M = conj(A) . W . A^T

    and M does not depend on the parameter either. So the two matrix products
    are done once here rather than once per parameter, and what is left per
    parameter is the elementwise sum against its own D -- turning 2*P*n^3 into
    2*n^3 + P*n^2, a factor of about 2n once there are more parameters than
    beams. It is the same finite sum in a different order, so the answer moves
    only by round-off.

    @param scale multiplies the amplitude of the output beam; folded into W so
      it costs nothing. Pass 1 where the formulation has no such factor.
    @param rv the already-scaled complex amplitude, for dI = 2 Re(conj(F) dF).
    */
    static void accumulate_grad_row(
      const cmat_t &A, const cmat_t &A_cjt, const cmat_t &G,
      size_t idx, complex_t scale,
      af::shared<cmat_t> const &Ds_kin,
      complex_t const &rv,
      mat_t &D_dyn)
    {
      const size_t n = A.accessor().n_columns();
      cmat_t W(af::mat_grid(n, n)),
        At(af::mat_grid(n, n)),   // A transposed
        cA(af::mat_grid(n, n));   // A conjugated, elementwise
      // G scaled by the output row of A and by A* first column: everything the
      // old per-parameter Hadamard product applied to V, gathered up before V
      // exists. The transposes are built in the same sweep, being free here.
      for (size_t i = 0; i < n; i++) {
        const complex_t a_i = A(idx, i)*scale;
        for (size_t j = 0; j < n; j++) {
          W(i, j) = a_i*G(i, j)*A_cjt(j, 0);
          At(j, i) = A(i, j);
          cA(i, j) = std::conj(A(i, j));
        }
      }
      // M = conj(A) . W . A^T, the two products that used to be A* D_p A once
      // per parameter. Done once, with no D in them at all.
      cmat_t M = af::matrix_multiply(
        cA.const_ref(),
        af::matrix_multiply(W.const_ref(), At.const_ref()).const_ref());
      for (size_t pi = 0; pi < Ds_kin.size(); pi++) {
        // what is left of the parameter's own contribution: M against D_p,
        // elementwise, which is dF/dp for this beam
        const complex_t *d = Ds_kin[pi].begin();
        const complex_t *m = M.begin();
        complex_t dp = 0;
        for (size_t k = 0; k < n*n; k++) {
          dp += m[k]*d[k];
        }
        // copy result to output (dI/dp -> |CI|^2)
        D_dyn(0, pi) = 2*(rv.real()*dp.real() + rv.imag()*dp.imag());
      }
    }

    af::shared<miller::index<> > indices;
    cmat_t A;
    mat3_t RMf;
    FloatType thickness;
  };


  enum {
    DYN_CALCULATOR_DEFAULT = DYN_MATRIX_DEFAULT,
    DYN_CALCULATOR_2013 = DYN_MATRIX_2013,
    DYN_CALCULATOR_2015 = DYN_MATRIX_2015
  };

  template <typename FloatType>
  class dyn_calculator_factory {
  public:
    ED_UTIL_TYPEDEFS;
    dyn_calculator_factory(int type);

    boost::shared_ptr<a_dyn_calculator<FloatType> > make(
      const af::shared<miller::index<> >& indices,
      const cmat_t& mat_Ug,
      const cart_t& K,
      const mat3_t& RMf,
      const cart_t& N,
      FloatType thickness) const;

    boost::shared_ptr<a_dyn_calculator<FloatType> > make(
      const af::shared<miller::index<> >& indices,
      const cart_t& K, const cart_t& N,
      FloatType thickness) const;
  private:
    int type;
  };

}}
