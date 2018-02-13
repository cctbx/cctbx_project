#ifndef CCTBX_MAPTBX_TARGET_AND_GRADIENTS_H
#define CCTBX_MAPTBX_TARGET_AND_GRADIENTS_H

#include <cctbx/maptbx/interpolation.h>

namespace cctbx {
  namespace maptbx {
    namespace target_and_gradients {
      namespace diffmap {

class compute {
public:
  compute(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<double, af::c_grid_padded<3> > const& map_target,
    af::const_ref<double, af::c_grid_padded<3> > const& map_current,
    double const& step,
    af::const_ref<scitbx::vec3<double> > const& sites_frac)
  {
    int nx = static_cast<int>(map_target.accessor().focus()[0]);
    int ny = static_cast<int>(map_target.accessor().focus()[1]);
    int nz = static_cast<int>(map_target.accessor().focus()[2]);
    af::versa<double, af::c_grid_padded<3> > diff_density_array;
    diff_density_array.resize(af::c_grid_padded<3>(nx,ny,nz), 0);
    int x1box=0   ;
    int x2box=nx-1;
    int y1box=0   ;
    int y2box=ny-1;
    int z1box=0   ;
    int z2box=nz-1;
    target_ = 0;
    double* diff_density_array_ = diff_density_array.begin();
    for(int kx = x1box; kx <= x2box; kx++) {
      int mx = scitbx::math::mod_positive(kx,nx);
      int mxny = mx*ny;
      for(int ky = y1box; ky <= y2box; ky++) {
        int my = scitbx::math::mod_positive(ky,ny);
        int mxnypmynz = (mxny+my)*nz;
        for(int kz = z1box; kz <= z2box; kz++) {
          int mz = scitbx::math::mod_positive(kz,nz);
          double diff = map_target[mxnypmynz+mz]-map_current[mxnypmynz+mz];
          target_ += (diff * diff); // (mx*NY+my)*NZ+mz
          diff_density_array_[mxnypmynz+mz] = -2 * diff;
        }
      }
    }
    // gradients
    af::const_ref<double, af::c_grid_padded<3> > diffd =
      diff_density_array.const_ref();
    cctbx::cartesian<> step_x = cctbx::cartesian<>(step,0,0);
    cctbx::cartesian<> step_y = cctbx::cartesian<>(0,step,0);
    cctbx::cartesian<> step_z = cctbx::cartesian<>(0,0,step);
    double two_step = 2*step;
    gradients_.resize(sites_frac.size(), scitbx::vec3<double>(0,0,0));
    for(std::size_t i_site=0;i_site<sites_frac.size();i_site++) {
      cctbx::fractional<> const& site_frac = sites_frac[i_site];
      cctbx::cartesian<> site_cart = unit_cell.orthogonalize(site_frac);
      cctbx::fractional<> sxp = unit_cell.fractionalize(site_cart+step_x);
      cctbx::fractional<> sxm = unit_cell.fractionalize(site_cart-step_x);
      cctbx::fractional<> syp = unit_cell.fractionalize(site_cart+step_y);
      cctbx::fractional<> sym = unit_cell.fractionalize(site_cart-step_y);
      cctbx::fractional<> szp = unit_cell.fractionalize(site_cart+step_z);
      cctbx::fractional<> szm = unit_cell.fractionalize(site_cart-step_z);
      double gx = (eight_point_interpolation(diffd, sxp) -
                   eight_point_interpolation(diffd, sxm)) / two_step;
      double gy = (eight_point_interpolation(diffd, syp) -
                   eight_point_interpolation(diffd, sym)) / two_step;
      double gz = (eight_point_interpolation(diffd, szp) -
                   eight_point_interpolation(diffd, szm)) / two_step;
      gradients_[i_site] = scitbx::vec3<double>(gx,gy,gz);
    }
  }

  double target() { return target_; }
  af::shared<scitbx::vec3<double> > gradients() { return gradients_; }

protected:
  double target_;
  af::shared<scitbx::vec3<double> > gradients_;
};

}}}} // namespace cctbx::maptbx::target_and_gradients::diffmap

namespace cctbx {
  namespace maptbx {
    namespace target_and_gradients {
      namespace simple {

template <typename FloatType=double>
class compute
{
public:

  compute(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<FloatType, af::c_grid_padded<3> > const& map_data,
    af::const_ref<scitbx::vec3<FloatType> > const& sites_cart,
    FloatType delta,
    af::const_ref<bool> const& selection)
  {
    gradients_.resize(sites_cart.size(), scitbx::vec3<FloatType>(0,0,0));
    target_ = 0;
    scitbx::vec3<FloatType>* res = gradients_.begin();
    for(std::size_t i_site=0;i_site<sites_cart.size();i_site++,res++) {
      if(selection[i_site]) {
        target_ += eight_point_interpolation(
          map_data,
          unit_cell.fractionalize(sites_cart[i_site]));
        scitbx::vec3<FloatType> piv = sites_cart[i_site];
        scitbx::vec3<FloatType> piv_d = piv;
        for(unsigned i_axis=0;i_axis<3;i_axis++) {
          FloatType densities[2];
          for(unsigned i_sign=0;i_sign<2;i_sign++) {
            piv_d[i_axis] = (i_sign == 0 ? piv[i_axis] + delta
                                         : piv[i_axis] - delta);
            fractional<FloatType> site_frac = unit_cell.fractionalize(piv_d);
            densities[i_sign] = eight_point_interpolation(
              map_data, site_frac);
          }
          piv_d[i_axis] = piv[i_axis];
          (*res)[i_axis] = (densities[0] - densities[1]) / (2 * delta);
        }
      }
    }
  }

  compute(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<FloatType, af::c_grid_padded<3> > const& map_data,
    af::const_ref<scitbx::vec3<FloatType> > const& sites_cart,
    af::const_ref<bool> const& selection,
    std::string const& interpolation)
  {
    gradients_.resize(sites_cart.size(), scitbx::vec3<FloatType>(0,0,0));
    af::c_grid_padded<3> a = map_data.accessor();
    scitbx::vec3<FloatType> step;
    for(unsigned i=0;i<3;i++) {
      step[i] = unit_cell.parameters()[i] / a.all()[i];
    }
    target_ = 0;
    for(std::size_t i_site=0;i_site<sites_cart.size();i_site++) {
      if(selection[i_site]) {
        af::tiny<FloatType, 4> result;
        if(interpolation == "linear") {
          result = eight_point_interpolation_with_gradients(
            map_data,
            unit_cell.fractionalize(sites_cart[i_site]),
            step);
        }
        else if(interpolation=="quadratic") {
          result = quadratic_interpolation_with_gradients(
            map_data,
            unit_cell.fractionalize(sites_cart[i_site]),
            step);
        }
        else if(interpolation=="tricubic") {
          result = tricubic_interpolation_with_gradients(
            map_data,
            unit_cell.fractionalize(sites_cart[i_site]),
            step);
        }
        else {
          throw std::runtime_error("Unknown interpolation mode.");
        }
        target_ += result[0];
        gradients_[i_site]=scitbx::vec3<FloatType>(
          result[1],result[2],result[3]);
      }
    }
  }

  FloatType target_;
  af::shared<scitbx::vec3<FloatType> > gradients_;

  FloatType target() { return target_; }
  af::shared<scitbx::vec3<FloatType> > gradients() { return gradients_; }
};


template <
  typename MapFloatType,
  typename SiteFloatType>
MapFloatType
target(
  uctbx::unit_cell const& unit_cell,
  af::const_ref<MapFloatType, af::c_grid_padded<3> > const& density_map,
  af::const_ref<scitbx::vec3<SiteFloatType> > const& sites_cart,
  af::const_ref<bool> const& selection)
{
  MapFloatType result = 0;
  for(std::size_t i_site=0;i_site<sites_cart.size();i_site++) {
    if(selection[i_site]) {
      result += eight_point_interpolation(
        density_map,
        unit_cell.fractionalize(sites_cart[i_site]));
    }
  }
  return result;
}

template <
  typename MapFloatType,
  typename SiteFloatType>
MapFloatType
target(
  uctbx::unit_cell const& unit_cell,
  af::const_ref<MapFloatType, af::c_grid_padded<3> > const& density_map,
  af::const_ref<scitbx::vec3<SiteFloatType> > const& sites_cart)
{
  MapFloatType result = 0;
  for(std::size_t i_site=0;i_site<sites_cart.size();i_site++) {
      result += eight_point_interpolation(
        density_map,
        unit_cell.fractionalize(sites_cart[i_site]));
  }
  return result;
}

template <
  typename MapFloatType,
  typename SiteFloatType>
MapFloatType
target(
  uctbx::unit_cell const& unit_cell,
  af::const_ref<MapFloatType, af::c_grid_padded<3> > const& density_map,
  af::const_ref<scitbx::vec3<SiteFloatType> > const& sites_cart,
  af::const_ref<std::size_t> const& selection)
{
  MapFloatType result = 0;
  for(std::size_t i_site=0;i_site<selection.size();i_site++) {
    result += eight_point_interpolation(
      density_map,
      unit_cell.fractionalize(sites_cart[selection[i_site]]));
  }
  return result;
}

template <
  typename MapFloatType,
  typename SiteFloatType>
MapFloatType
target(
  af::const_ref<MapFloatType, af::c_grid<3> > const& density_map,
  af::const_ref<scitbx::vec3<SiteFloatType> > const& sites_frac)
{
  MapFloatType result = 0;
  for(std::size_t i_site=0;i_site<sites_frac.size();i_site++) {
    result += eight_point_interpolation(
      density_map,
      sites_frac[i_site]);
  }
  return result;
}

template <
  typename MapFloatType,
  typename SiteFloatType>
af::shared<MapFloatType>
target_per_site(
  uctbx::unit_cell const& unit_cell,
  af::const_ref<MapFloatType, af::c_grid_padded<3> > const& density_map,
  af::const_ref<scitbx::vec3<SiteFloatType> > const& sites_cart)
{
  af::shared<MapFloatType> result(
    sites_cart.size(), af::init_functor_null<MapFloatType>());
  for(std::size_t i_site=0;i_site<sites_cart.size();i_site++) {
    result[i_site] = eight_point_interpolation(
      density_map,
      unit_cell.fractionalize(sites_cart[i_site]));
  }
  return result;
}

template <
  typename MapFloatType,
  typename SiteFloatType>
af::shared<scitbx::vec3<SiteFloatType> >
gradients(
  uctbx::unit_cell const& unit_cell,
  af::const_ref<MapFloatType, af::c_grid_padded<3> > const& density_map,
  af::const_ref<scitbx::vec3<SiteFloatType> > const& sites_cart,
  SiteFloatType delta,
  af::const_ref<bool> const& selection)
{
  CCTBX_ASSERT(delta > 0);
  typedef scitbx::vec3<SiteFloatType> v3_t;
  af::shared<v3_t> result(sites_cart.size(), af::init_functor_null<v3_t>());
  v3_t* res = result.begin();
  for(std::size_t i_site=0;i_site<sites_cart.size();i_site++,res++) {
    result[i_site]=scitbx::vec3<SiteFloatType>(0,0,0);
    if(selection[i_site]) {
      v3_t piv = sites_cart[i_site];
      v3_t piv_d = piv;
      for(unsigned i_axis=0;i_axis<3;i_axis++) {
        MapFloatType densities[2];
        for(unsigned i_sign=0;i_sign<2;i_sign++) {
          piv_d[i_axis] = (i_sign == 0 ? piv[i_axis] + delta
                                       : piv[i_axis] - delta);
          fractional<SiteFloatType> site_frac = unit_cell.fractionalize(piv_d);
          densities[i_sign] = eight_point_interpolation(
            density_map, site_frac);
        }
        piv_d[i_axis] = piv[i_axis];
        (*res)[i_axis] = (densities[0] - densities[1]) / (2 * delta);
      }
    }
  }
  return result;
}

}}}} // namespace cctbx::maptbx::target_and_gradients::simple

#endif // GUARD
