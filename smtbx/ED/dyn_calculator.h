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
      FloatType thickness)
      : indices(indices),
      K(K),
      thickness(thickness)
    {}

    virtual ~a_dyn_calculator() {}

    // mat_Ug will be NOT be affected - deep copied
    a_dyn_calculator& reset(const cmat_t& m, const mat3_t &RMf, const cart_t &N) {
      A = m.deep_copy();
      this->RMf = RMf;
      this->N = N;
      return build();
    }

    // mat_Ug will be NOT be affected - deep copied
    a_dyn_calculator& reset(const af::shared<miller::index<> > &indices_,
      const cmat_t& m, const mat3_t& RMf, const cart_t& N) {
      A = m.deep_copy();
      indices = indices_;
      this->RMf = RMf;
      this->N = N;
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
    const cart_t K;
  protected:
    af::shared<miller::index<> > indices;
    cmat_t A;
    mat3_t RMf;
    cart_t N;
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
      const cart_t& K,
      FloatType thickness) const;
  private:
    int type;
  };
  
}}
