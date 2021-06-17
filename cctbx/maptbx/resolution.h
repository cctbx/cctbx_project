#ifndef RESOLUTION_H
#define RESOLUTION_H

#include <cctbx/miller.h>

namespace cctbx { namespace maptbx {

//! Thin box worth of hkl starting from highest d_min with some step. For each
//! step calculate CC between original map and map computed using truncated box.
class d99 {
public:
  double d_cut_;
  d99(
    af::const_ref<std::complex<double> > const& f,
    af::const_ref<double> const& d_spacings,
    af::const_ref<miller::index<> > const& hkl,
    double const& cutoff)
  {
    double d_min =  1.e+9;
    double d_max = -1.e+9;
    for(int i = 0; i < d_spacings.size(); i++) {
      double d = d_spacings[i];
      if(d < d_min) d_min = d;
      if(d > d_max) d_max = d;
    }
    CCTBX_ASSERT(f.size() == d_spacings.size());
    CCTBX_ASSERT(f.size() == hkl.size());
    af::shared<double> data_sq_w(f.size());
    // Precompute intermediate arrays
    double sum_all = 0.;
    for(int i = 0; i < f.size(); i++) {
      miller::index<> h = hkl[i];
      double w;
      if(h[2]==0) { w = 1.; }
      else        { w = 2.; }
      double f_abs = std::abs(f[i]);
      double f_abs_sq = f_abs*f_abs;
      double w_f_abs_sq = w*f_abs_sq;
      sum_all += w_f_abs_sq;
      data_sq_w[i] = w_f_abs_sq;
    }
    CCTBX_ASSERT(sum_all != 0.);
    //
    double step=0.1;
    double d=d_min;
    d_cut_ = -1;
    while(d<d_max) {
      d+=step;
      double sum_num = 0.;
      for(int i = 0; i < f.size(); i++) {
        if(d_spacings[i]>d) {
          sum_num += data_sq_w[i];
        }
      }
      double cc = std::sqrt(sum_num/sum_all);
      if(cc < cutoff) {
        d_cut_ = d;
        break;
      }
    }
    //
    step=0.01;
    d=d_cut_-0.2;
    while(d<d_max) {
      d+=step;
      double sum_num = 0.;
      for(int i = 0; i < f.size(); i++) {
        if(d_spacings[i]>d) {
          sum_num += data_sq_w[i];
        }
      }
      double cc = std::sqrt(sum_num/sum_all);
      if(cc < cutoff) {
        d_cut_ = d;
        break;
      }
    }
    //
  }

  double d_min() {return d_cut_;}

};

}} // namespace cctbx::maptbx

#endif
