#ifndef FT_ANALYTICAL_1D_POINT_SCATTERER_H
#define FT_ANALYTICAL_1D_POINT_SCATTERER_H

#include <scitbx/array_family/accessors/c_grid.h>
#include <cctbx/xray/sampling_base.h>

namespace cctbx { namespace maptbx {

template <typename FloatType>
class ft_analytical_1d_point_scatterer_at_origin {
public:
  af::shared<FloatType> distances_;
  af::shared<FloatType> rho_;
  af::shared<FloatType> cos_table;
  int N;
  FloatType st;

  ft_analytical_1d_point_scatterer_at_origin(int const& N_)
  :
  N(N_)
  {
    st = 2*scitbx::constants::pi/N;
    for(int i=0; i<N; i++) cos_table.push_back(std::cos(i*st));
  }

  void compute(
    af::shared<miller::index<> > miller_indices,
    FloatType const& step,
    FloatType const& left,
    FloatType const& right,
    af::shared<FloatType> const& u_frac)
  {
    distances_ = af::shared<FloatType>();
    rho_ = af::shared<FloatType>();
    FloatType two_pi = 2*scitbx::constants::pi;
    FloatType r = left;
    while(r<=right) {
      FloatType two_pi_r = two_pi * r;
      FloatType rho__ = 0;
      FloatType x=u_frac[0];
      FloatType y=u_frac[1];
      FloatType z=u_frac[2];
      for(std::size_t i=0; i<miller_indices.size(); i++) {
        miller::index<> mi = miller_indices[i];
        FloatType arg = two_pi_r*(x*mi[0]+y*mi[1]+z*mi[2]);
        // use standard cosine
        //rho__+=std::cos(arg);
        // use table
        if(arg<0) arg = std::abs(arg);
        if(arg>two_pi) arg = arg - int(arg/two_pi)*two_pi;
        arg = arg/st;
        // this is faster but may be dangerous (is it really dangerous?)
        //int k = int(arg);
        //int k1= k+1;
        int k = scitbx::math::mod_positive(int(arg),N);
        FloatType y=cos_table[k];
        rho__ += (y+(cos_table[scitbx::math::mod_positive(k+1,N)]-y)*(arg-k));
        //
      }
      distances_.push_back(r);
      rho_.push_back(rho__);
      r+=step;
    }
  }

  af::shared<FloatType> rho() { return rho_; }

  af::shared<FloatType> distances() { return distances_; }

};

}} // namespace cctbx::maptbx

#endif // FT_ANALYTICAL_1D_POINT_SCATTERER_H
