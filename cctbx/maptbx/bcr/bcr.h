#ifndef BCR_H
#define BCR_H

#include <boost/python/list.hpp>
#include <cctbx/maptbx/interpolation.h> // indirect import
#include <cctbx/xray/scatterer.h>
#include <cctbx/adptbx.h>
#include <cctbx/xray/sampling_base.h>
#include <iotbx/pdb/hierarchy.h>

namespace cctbx { namespace maptbx {

using scitbx::constants::pi;
using scitbx::constants::four_pi_sq;
using scitbx::constants::four_pi;

template <typename FloatType=double,
          typename XrayScattererType=cctbx::xray::scatterer<> >
class bcr_scatterer
{
  public:
    XrayScattererType const& scatterer;
    FloatType radius;
    FloatType resolution;
    af::shared<FloatType> mu;
    af::shared<FloatType> kappa;
    af::shared<FloatType> nu;
    af::shared<FloatType> musq;
    af::shared<FloatType> kappi;

    bcr_scatterer(
      XrayScattererType const& scatterer_,
      FloatType radius_,
      FloatType resolution_,
      af::shared<FloatType> mu_,
      af::shared<FloatType> kappa_,
      af::shared<FloatType> nu_,
      af::shared<FloatType> musq_,
      af::shared<FloatType> kappi_
      )
    :
      radius(radius_),
      resolution(resolution_),
      mu(mu_),
      kappa(kappa_),
      nu(nu_),
      musq(musq_),
      kappi(kappi_),
      scatterer(scatterer_)
    {
      CCTBX_ASSERT(mu.size()==kappa.size());
      CCTBX_ASSERT(mu.size()==nu.size());
      CCTBX_ASSERT(mu.size()==musq.size());
      CCTBX_ASSERT(mu.size()==kappi.size());
    }

    XrayScattererType const& get_scatterer() const { return scatterer; }

};

template <typename FloatType=double,
          typename XrayScattererType=cctbx::xray::scatterer<> >
class bcr_model
{
  public:
    XrayScattererType scatterer;
    af::shared<FloatType> B;
    af::shared<FloatType> C;
    af::shared<FloatType> R;
    FloatType b_iso;

    bcr_model() {
      scatterer = cctbx::xray::scatterer<> ();
      B.fill(0.0);
      C.fill(0.0);
      R.fill(0.0);
      b_iso = cctbx::adptbx::u_as_b(scatterer.u_iso);
    }

    bcr_model(
      XrayScattererType     const& scatterer_,
      af::shared<FloatType> const& B_,
      af::shared<FloatType> const& C_,
      af::shared<FloatType> const& R_
      )
    :
      scatterer(scatterer_),
      B(B_),
      C(C_),
      R(R_),
      b_iso(cctbx::adptbx::u_as_b(scatterer.u_iso))
    {}
};

template <typename FloatType=double,
          typename XrayScattererType=cctbx::xray::scatterer<> >
class calculator
{
  public:
    bcr_model<FloatType> bcrm;

  // Default constructor
  calculator() {}

  // Exp table provided and will be used (available only in C++)
  calculator(
    bcr_model<FloatType> const& bcrm_,
    cctbx::xray::detail::exponent_table<FloatType>& exp_table_)
  :
    bcrm(bcrm_),
    exp_table(&exp_table_),
    use_exp_table(true)
  {}

  // std::exp will be used. Available in C++ and Python
  calculator(
    bcr_model<FloatType> const& bcrm_)
  :
    bcrm(bcrm_),
    use_exp_table(false)
  {}

  FloatType atom_radius() const
  {
    return 3.;
  }

  FloatType rho(FloatType const& r) const
  {
    FloatType result=0;
    for(std::size_t i=0; i<bcrm.B.size(); i++) {
      FloatType B_i = bcrm.B[i];
      FloatType C_i = bcrm.C[i];
      FloatType R_i = bcrm.R[i];
      FloatType B_plus = B_i + bcrm.b_iso;
      FloatType result_g=0;
      FloatType result_c=0;
      if(std::abs(bcrm.R[i])<1.e-6) { // gauss
        FloatType fpob = 4 * pi / B_plus;
        FloatType XXX1 = -fpob * pi * r*r;//---
        //std::cout<<"XXX1 "<<XXX1<<" fpob "<<fpob<<"r "<<r<<std::endl;
        result_g += C_i * std::pow(fpob,1.5) * myexp(XXX1);
      } //  end gauss
      else { // chi
        FloatType mfpsob = (-four_pi_sq)/B_plus;
        if(std::abs(r) < 1.e-6) {
          FloatType XXX2 = mfpsob*R_i*R_i;//---
          //std::cout<<"XXX2 "<<XXX2<<std::endl;
          result_c += std::pow(four_pi/B_plus,1.5) * myexp(XXX2);
        }
        else {
          FloatType XXX3 = mfpsob*std::pow(r-R_i,2);
          FloatType XXX4 = mfpsob*std::pow(r+R_i,2);
          //std::cout<<"XXX3 "<<XXX3<<std::endl;
          //std::cout<<"XXX4 "<<XXX4<<std::endl;
          FloatType e1 = myexp(XXX3);
          FloatType e2 = myexp(XXX4);
          result_c += std::pow(four_pi*B_plus,-0.5) * (1/(r*R_i)) * (e1-e2);
        }
        result_c *= C_i;
      } // end chi
      result += (result_g+result_c);
    }
    return result;
  }

  protected:
    cctbx::xray::detail::exponent_table<FloatType>* exp_table;
  private:
    bool use_exp_table;
    inline FloatType myexp(FloatType const& x) const
    {
      //std::cout<<x<<std::endl;
      if(use_exp_table) return (*exp_table)(x);
      else              return std::exp(x);
    }
};

template <typename FloatType=double,
          typename XrayScattererType=cctbx::xray::scatterer<> >
class image {
public:
  af::shared<double> cc_, d_, d_inv_;
  image(
    cctbx::uctbx::unit_cell const& unit_cell,
    boost::python::list const& bcr_models,
    int const& step)
  {
    int il=0;
    int ir=step;
    for(std::size_t i=0; i<boost::python::len(bcr_models); i++) {
       bcr_model<FloatType> bcrm =
         boost::python::extract<bcr_model<FloatType> >(bcr_models[i])();
       XrayScattererType scatteter = bcrm.scatterer;

       std::cout<<"HERE"<<std::endl;
       //std::cout<<bcrm.rho(1.3)<<std::endl;

       calculator<FloatType> calc1 = calculator<FloatType>(bcrm);
       std::cout<<calc1.rho(1.3)<<std::endl;

       cctbx::xray::detail::exponent_table<FloatType> exp_table(-100);
       calculator<FloatType> calc2 = calculator<FloatType>(bcrm, exp_table);
       std::cout<<calc2.rho(1.3)<<std::endl;

    }

  }

  af::shared<double> cc()     {return cc_;}
  af::shared<double> d()      {return d_;}
  af::shared<double> d_inv()  {return d_inv_;}

};


}} // namespace cctbx::maptbx

#endif
