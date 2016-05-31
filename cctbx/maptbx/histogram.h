#ifndef CCTBX_MAPTBX_HISTOGRAM_H
#define CCTBX_MAPTBX_HISTOGRAM_H

#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/math/utils.h>

namespace cctbx { namespace maptbx {

//! Compute map histogram and cumulative functions
class histogram {
public:
  af::shared<double> values_;
  af::shared<double> c_values_;
  af::shared<double> v_values_;
  af::shared<double> arguments_;
  double bin_width_;
  int n_bins_;
  histogram(
    af::const_ref<double, af::c_grid<3> > const& map,
    int const& n_bins__,
    double min=-1,
    double max=-1)
  :
  n_bins_(n_bins__)
  {
    int nx = map.accessor()[0];
    int ny = map.accessor()[1];
    int nz = map.accessor()[2];
    if(std::abs(std::abs(min)-std::abs(max))<1.e-9) {
      max = af::max(map);
      min = af::min(map);
    }
    double size = map.size();
    CCTBX_ASSERT(size > 0);
    CCTBX_ASSERT(n_bins_ > 0);
    values_.resize(n_bins_, 0);
    arguments_.resize(n_bins_, 0);
    bin_width_ = (max-min)/n_bins_;
    CCTBX_ASSERT(bin_width_ > 0);
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
          double rho = map(i,j,k);
          int index = scitbx::math::nearest_integer((rho-min)/bin_width_);
          if(index<0) index=0;
          if(index>=n_bins_) index=n_bins_-1;
          values_[index] += 1;
    }}}
    for (int i = 0; i < n_bins_; i++) {
      values_[i] = values_[i] / size;
    }
    c_values_.resize(n_bins_, 0);
    for (int i = 0; i < n_bins_; i++) {
      double sum = 0;
      for (int j = 0; j < i; j++) {
        sum += values_[j];
      }
      c_values_[i] = sum;
    }
    v_values_.resize(n_bins_, 0);
    for (int i = 0; i < n_bins_; i++) {
      v_values_[i] = 1. - c_values_[i];
      arguments_[i] = min+i*bin_width_;
    }
  }

  histogram(
    af::const_ref<double> const& map,
    int const& n_bins__,
    double min=-1,
    double max=-1)
  :
  n_bins_(n_bins__)
  {
    if(std::abs(std::abs(min)-std::abs(max))<1.e-9) {
      max = af::max(map);
      min = af::min(map);
    }
    double size = map.size();
    CCTBX_ASSERT(size > 0);
    CCTBX_ASSERT(n_bins_ > 0);
    values_.resize(n_bins_, 0);
    arguments_.resize(n_bins_, 0);
    bin_width_ = (max-min)/n_bins_;
    CCTBX_ASSERT(bin_width_ > 0);
    for (int i = 0; i < map.size(); i++) {
      double rho = map[i];
      int index = scitbx::math::nearest_integer((rho-min)/bin_width_);
      if(index<0) index=0;
      if(index>=n_bins_) index=n_bins_-1;
      values_[index] += 1;
    }
    for (int i = 0; i < n_bins_; i++) {
      values_[i] = values_[i] / size;
    }
    c_values_.resize(n_bins_, 0);
    for (int i = 0; i < n_bins_; i++) {
      double sum = 0;
      for (int j = 0; j < i; j++) {
        sum += values_[j];
      }
      c_values_[i] = sum;
    }
    v_values_.resize(n_bins_, 0);
    for (int i = 0; i < n_bins_; i++) {
      v_values_[i] = 1. - c_values_[i];
      arguments_[i] = min+i*bin_width_;
    }
  }

  af::shared<double> values()    {return values_;}
  af::shared<double> c_values()  {return c_values_;}
  af::shared<double> v_values()  {return v_values_;}
  double bin_width()             {return bin_width_;}
  af::shared<double> arguments() {return arguments_;}
  int n_bins()                   {return n_bins_;}

};

}} // namespace cctbx::maptbx

#endif
