#ifndef FSC_H
#define FSC_H

#include <cctbx/miller.h>
#include <cctbx/maptbx/utils.h>

namespace cctbx { namespace maptbx {

class fsc {
public:
  af::shared<double> cc_, d_, d_inv_;
  fsc(
    af::const_ref<std::complex<double> > const& f1,
    af::const_ref<std::complex<double> > const& f2,
    af::const_ref<double> const& d_spacings,
    int const& step)
  {
    CCTBX_ASSERT(f1.size() == d_spacings.size());
    CCTBX_ASSERT(f1.size() == f2.size());
    int il=0;
    int ir=step;
    int size = f1.size();
    while(ir<size) {
      af::shared<std::complex<double> > f1_i;
      af::shared<std::complex<double> > f2_i;
      double d_mean = 0;
      for(int i = il; i < ir; i++) {
        f1_i.push_back(f1[i]);
        f2_i.push_back(f2[i]);
        d_mean += d_spacings[i];
      }
      d_mean = d_mean/step;
      d_.push_back(d_mean);
      d_inv_.push_back(1./d_mean);
      cc_.push_back( cc_complex_complex<std::complex<double>, double>(
        f1_i.const_ref(), f2_i.const_ref()) );
      il += step;
      ir += step;
      if(ir>=size) break;
    }
  }

  af::shared<double> cc()     {return cc_;}
  af::shared<double> d()      {return d_;}
  af::shared<double> d_inv()  {return d_inv_;}

};

}} // namespace cctbx::maptbx

#endif
