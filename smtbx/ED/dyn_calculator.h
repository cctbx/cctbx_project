#pragma once
#include <smtbx/ED/ed.h>
#include <boost/shared_ptr.hpp>
namespace smtbx { namespace ED
{
  using namespace cctbx;

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
      for (size_t i = 0; i < n; i++) {
        const complex_t a_i = A(idx, i)*scale;
        for (size_t j = 0; j < n; j++) {
          W(i, j) = a_i*G(i, j)*A_cjt(j, 0);
          At(j, i) = A(i, j);
          cA(i, j) = std::conj(A(i, j));
        }
      }
      cmat_t M = af::matrix_multiply(
        cA.const_ref(),
        af::matrix_multiply(W.const_ref(), At.const_ref()).const_ref());
      for (size_t pi = 0; pi < Ds_kin.size(); pi++) {
        const complex_t *d = Ds_kin[pi].begin();
        const complex_t *m = M.begin();
        complex_t dp = 0;
        for (size_t k = 0; k < n*n; k++) {
          dp += m[k]*d[k];
        }
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
