#ifndef RESOLUTION_H
#define RESOLUTION_H

#include <cctbx/miller.h>

namespace cctbx { namespace maptbx {

//! Thin box worth of hkl starting from highest d_min with some step. For each
//! step calculate CC between original map and map computed using truncated box.
class d99 {
public:
  double d_min_cc9_, d_min_cc99_, d_min_cc999_, d_min_cc9999_, d_min_cc99999_;
  double d_min_cc999999_;
  d99(
    af::const_ref<std::complex<double> > const& f,
    af::const_ref<double> const& d_spacings,
    af::const_ref<miller::index<> > const& hkl,
    double const& d_min,
    double const& d_max)
  {
    CCTBX_ASSERT(f.size() == d_spacings.size());
    CCTBX_ASSERT(f.size() == hkl.size());
    d_min_cc9_=-1., d_min_cc99_=-1., d_min_cc999_=-1., d_min_cc9999_=-1;
    d_min_cc99999_=-1, d_min_cc999999_=-1.;
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
      sum_all += (w*f_abs_sq);
      data_sq_w[i] = w*f_abs_sq;
    }
    CCTBX_ASSERT(sum_all != 0.);
    //
    double step=1.0;
    double d=d_min;
    double d_cut=0;
    double tmp0=0;
    while(d<d_max) {
      d+=step;
      double sum_num = 0.;
      for(int i = 0; i < f.size(); i++) {
        if(d_spacings[i]>d) {
          sum_num += data_sq_w[i];
        }
      }
      double cc = std::sqrt(sum_num/sum_all);
      if(cc<0.99-0.01) {
        d_cut = d;
        tmp0 = cc;
        break;
      }
    }
    //
    double diff=999.;
    int i_cut=0;
    double tmp = 0;
    for(int i = 0; i < f.size(); i++) {
      double diff_ = std::abs(d_spacings[i]-d_cut);
      if(diff_<diff) {
        diff = diff_;
        i_cut = i;
        tmp=d_spacings[i];
      }
    }
    //
    double sum_num_const=0;
    for(int i = i_cut; i < f.size(); i++) {
      sum_num_const += data_sq_w[i];
    }
    //
    int inc = 100;
    int il=i_cut-inc;
    int ir=i_cut;
    double sum_num_roll = sum_num_const;
    while(il>0) {
      double sum_num_range = 0.;
      for(int i = il; i < ir; i++) {
        sum_num_range += data_sq_w[i];
      }
      sum_num_roll += sum_num_range;
      double cc = std::sqrt(sum_num_roll/sum_all);
      d = d_spacings[il];
      if(d_min_cc999999_<0. && cc>0.999999) d_min_cc999999_=d;
      if(d_min_cc99999_<0. && cc>0.99999) d_min_cc99999_=d;
      if(d_min_cc9999_<0. && cc>0.9999) d_min_cc9999_=d;
      if(d_min_cc999_<0. && cc>0.999) d_min_cc999_=d;
      if(d_min_cc99_<0. && cc>0.99) d_min_cc99_=d;
      il = il-inc;
      ir = ir-inc;
      if(cc>0.999999) break;
    }
  }

  double d_min_cc9()      {return d_min_cc9_;}
  double d_min_cc99()     {return d_min_cc99_;}
  double d_min_cc999()    {return d_min_cc999_;}
  double d_min_cc9999()   {return d_min_cc9999_;}
  double d_min_cc99999()  {return d_min_cc99999_;}
  double d_min_cc999999() {return d_min_cc999999_;}

};

}} // namespace cctbx::maptbx

#endif
