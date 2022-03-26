#ifndef BCR_H
#define BCR_H

#include <boost/python/list.hpp>
#include <cctbx/maptbx/interpolation.h> // indirect import

namespace cctbx { namespace maptbx {

template <typename FloatType=double>
class bcr_model
{
  public:
    scitbx::vec3<FloatType> site_frac;
    af::shared<FloatType> b;
    af::shared<FloatType> c;
    af::shared<FloatType> r;

    bcr_model() {
      site_frac.fill(0.0);
      b.fill(0.0);
      c.fill(0.0);
      r.fill(0.0);
    }

    bcr_model(
      scitbx::vec3<FloatType> const& site_frac_,
      af::shared<FloatType> const& b_,
      af::shared<FloatType> const& c_,
      af::shared<FloatType> const& r_
      )
    :
      site_frac(site_frac_),
      b(b_),
      c(c_),
      r(r_)
    {}
};

template <typename FloatType=double>
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
       bcr_model<double> bcrm =
         boost::python::extract<bcr_model<double> >(bcr_models[i])();
       cctbx::fractional<double> site_frac = bcrm.site_frac;
    }

  }

  af::shared<double> cc()     {return cc_;}
  af::shared<double> d()      {return d_;}
  af::shared<double> d_inv()  {return d_inv_;}

};


}} // namespace cctbx::maptbx

#endif
