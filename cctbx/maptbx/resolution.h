#ifndef RESOLUTION_H
#define RESOLUTION_H

#include <cctbx/miller.h>

namespace cctbx { namespace maptbx {

//! Thin box worth of hkl starting from highest d_min with some step. For each
//! step calculate CC between original map and map computed using truncated box.
class d99 {
public:
  af::shared<double> ccs_;
  af::shared<double> d_mins_;
  double d_min_cc9_, d_min_cc99_, d_min_cc999_;
  d99(
    af::const_ref<std::complex<double> > const& f,
    af::const_ref<double> const& d_spacings,
    af::const_ref<miller::index<> > const& hkl,
    double const& d_min,
    double const& d_max)
  {
    CCTBX_ASSERT(f.size() == d_spacings.size());
    CCTBX_ASSERT(f.size() == hkl.size());
    d_min_cc9_=-1., d_min_cc99_=-1., d_min_cc999_=-1.;
    af::shared<double> w(f.size());
    af::shared<double> data_sq_w(f.size());
    double sum_all = 0.;
    for(int i = 0; i < f.size(); i++) {
      miller::index<> h = hkl[i];
      if(h[2]==0) { w[i] = 1.; }
      else        { w[i] = 2.; }
      double f_abs = std::abs(f[i]);
      double f_abs_sq = f_abs*f_abs;
      sum_all += (w[i]*f_abs_sq);
      data_sq_w[i] = w[i]*f_abs_sq;
    }
    CCTBX_ASSERT(sum_all != 0.);
    double step=0.01;
    double d=d_min;
    while(d<d_max) {
      d+=step;
      double sum_num = 0.;
      for(int i = 0; i < f.size(); i++) {
        if(d_spacings[i]>d) {
          sum_num += data_sq_w[i];
        }
      }
      double cc = std::sqrt(sum_num/sum_all);
      ccs_.push_back(cc);
      d_mins_.push_back(d);
      if(d_min_cc999_<0. && cc<0.999) d_min_cc999_=d;
      if(d_min_cc99_<0. && cc<0.99) {
        d_min_cc99_=d;
        step=0.1;
      }
      if(d_min_cc9_<0 && cc<0.9) d_min_cc9_=d;
      if(cc<0.9) break; // Perhaps this should be a paraameter.
    }
  }

  af::shared<double> ccs()    {return ccs_;}
  af::shared<double> d_mins() {return d_mins_;}
  double d_min_cc9()   {return d_min_cc9_;}
  double d_min_cc99()  {return d_min_cc99_;}
  double d_min_cc999() {return d_min_cc999_;}

};

}} // namespace cctbx::maptbx

#endif
