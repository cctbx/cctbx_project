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

using scitbx::mat3;

//- Magnification begin --------------------------------------------------------

template <
  typename MapFloatType,
  typename SiteFloatType>
MapFloatType
magnification_isotropic(
  uctbx::unit_cell const& unit_cell,
  af::const_ref<MapFloatType, af::c_grid_padded<3> > const& density_map,
  af::const_ref<scitbx::vec3<SiteFloatType> > const& sites_cart)
{
  MapFloatType t_best = 0;
  for(std::size_t i_site=0;i_site<sites_cart.size();i_site++) {
    t_best += tricubic_interpolation(
      density_map,
      unit_cell.fractionalize(sites_cart[i_site]));
  }
  MapFloatType m_min  = 0.9;
  MapFloatType m_max  = 1.1;
  MapFloatType inc    = 0.0001;
  MapFloatType m_best = 1.0;
  while(m_min<=m_max) {
    MapFloatType t = 0;
    for(std::size_t i_site=0;i_site<sites_cart.size();i_site++) {
      t += eight_point_interpolation( //tricubic_interpolation(
        density_map,
        unit_cell.fractionalize(sites_cart[i_site])*m_min);
    }
    if(t>t_best) {
      t_best = t;
      m_best = m_min;
    }
    m_min += inc;
  }
  return m_best;
}

// Too slow to be practical in most cases. Also, inc needs to be 0.0001 which
// will make it impossible.
template <
  typename MapFloatType,
  typename SiteFloatType>
scitbx::vec3<MapFloatType>
magnification_anisotropic(
  uctbx::unit_cell const& unit_cell,
  af::const_ref<MapFloatType, af::c_grid_padded<3> > const& density_map,
  af::const_ref<scitbx::vec3<SiteFloatType> > const& sites_cart)
{
  MapFloatType t_best = 0;
  for(std::size_t i_site=0;i_site<sites_cart.size();i_site++) {
    t_best += tricubic_interpolation(
      density_map,
      unit_cell.fractionalize(sites_cart[i_site]));
  }
  MapFloatType m_min  = 0.9;
  MapFloatType m_max  = 1.1;
  MapFloatType inc    = 0.01;
  MapFloatType m1_best = 1.0;
  MapFloatType m2_best = 1.0;
  MapFloatType m3_best = 1.0;
  MapFloatType m1 = m_min;
  while(m1<=m_max) {
    MapFloatType m2 = m_min;
    while(m2<=m_max) {
      MapFloatType m3 = m_min;
      while(m3<=m_max) {
        MapFloatType t = 0;
        for(std::size_t i_site=0;i_site<sites_cart.size();i_site++) {
          scitbx::vec3<SiteFloatType> s = sites_cart[i_site];
          scitbx::vec3<SiteFloatType> sm =
            scitbx::vec3<SiteFloatType>(s[0]*m1,s[1]*m2,s[2]*m3);
          t += eight_point_interpolation(
            density_map,
            unit_cell.fractionalize(sm));
        }
        if(t>t_best) {
          t_best = t;
          m1_best = m1;
          m2_best = m2;
          m3_best = m3;
        }
        m3 += inc;
      }
      m2 += inc;
    }
    m1 += inc;
  }
  return scitbx::vec3<MapFloatType>(m1_best, m2_best, m3_best);
}

template <typename FloatType=double>
class magnification
{
public:

  magnification(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<FloatType, af::c_grid_padded<3> > const& map_data,
    af::const_ref<scitbx::vec3<FloatType> > const& sites_cart,
    mat3<FloatType> const& K)
  {
    gradients_.resize(9, 0);
    af::c_grid_padded<3> a = map_data.accessor();
    scitbx::vec3<FloatType> step;
    for(unsigned i=0;i<3;i++) {
      step[i] = unit_cell.parameters()[i] / a.all()[i];
    }
    target_ = 0;
    for(std::size_t i_site=0;i_site<sites_cart.size();i_site++) {
      scitbx::vec3<FloatType> site_cart = sites_cart[i_site];
      FloatType x = site_cart[0];
      FloatType y = site_cart[1];
      FloatType z = site_cart[2];
      af::tiny<FloatType, 4> result = tricubic_interpolation_with_gradients(
        map_data,
        unit_cell.fractionalize(K * site_cart),
        step);
      target_ += result[0];
      scitbx::vec3<FloatType> gXYZ =
        unit_cell.orthogonalize(scitbx::vec3<FloatType>(
          result[1],result[2],result[3]));
      FloatType gx = gXYZ[0];
      FloatType gy = gXYZ[1];
      FloatType gz = gXYZ[2];

      gradients_[0] += gx*x; // gk11
      gradients_[1] += gx*y; // gk12
      gradients_[2] += gx*z; // gk13
      gradients_[3] += gy*x; // gk21
      gradients_[4] += gy*y; // gk22
      gradients_[5] += gy*z; // gk23
      gradients_[6] += gz*x; // gk31
      gradients_[7] += gz*y; // gk32
      gradients_[8] += gz*z; // gk33
    }
  }

  magnification(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<FloatType, af::c_grid_padded<3> > const& map_data,
    af::const_ref<scitbx::vec3<FloatType> > const& sites_cart,
    scitbx::vec3<FloatType> const& K)
  {
    gradients_.resize(3, 0);
    af::c_grid_padded<3> a = map_data.accessor();
    scitbx::vec3<FloatType> step;
    for(unsigned i=0;i<3;i++) {
      step[i] = unit_cell.parameters()[i] / a.all()[i];
    }
    target_ = 0;
    for(std::size_t i_site=0;i_site<sites_cart.size();i_site++) {
      scitbx::vec3<FloatType> s = sites_cart[i_site];
      scitbx::vec3<FloatType> sm =
        scitbx::vec3<FloatType>(s[0]*K[0],s[1]*K[1],s[2]*K[2]);
      af::tiny<FloatType, 4> result = tricubic_interpolation_with_gradients(
        map_data,
        unit_cell.fractionalize(sm),
        step);
      target_ += result[0];
      scitbx::vec3<FloatType> gXYZ =
        unit_cell.orthogonalize(scitbx::vec3<FloatType>(
          result[1],result[2],result[3]));
      gradients_[0] += gXYZ[0]*s[0]; // gk11
      gradients_[1] += gXYZ[1]*s[1]; // gk22
      gradients_[2] += gXYZ[2]*s[2]; // gk33
    }
  }

  FloatType target_;
  af::shared<FloatType> gradients_;

  FloatType target() { return target_; }
  af::shared<FloatType> gradients() { return gradients_; }
};

//- Magnification end ----------------------------------------------------------

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

// Binary score
template <
  typename MapFloatType,
  typename SiteFloatType>
af::tiny<SiteFloatType, 2>
score(
  uctbx::unit_cell const& unit_cell,
  af::const_ref<MapFloatType, af::c_grid_padded<3> > const& density_map,
  af::const_ref<scitbx::vec3<SiteFloatType> > const& sites_cart,
  af::const_ref<std::size_t> const& selection,
  af::shared<af::tiny<std::size_t, 2> > bonded_pairs,
  af::const_ref<MapFloatType > const& weights,
  af::const_ref<MapFloatType > const& max_map_value,
  af::const_ref<MapFloatType > const& min_map_value)
{
  MapFloatType status = 1;
  MapFloatType result = 0;
  for(std::size_t i_site=0;i_site<selection.size();i_site++) {
    std::size_t i = selection[i_site];
    MapFloatType mv = eight_point_interpolation(
      density_map,
      unit_cell.fractionalize(sites_cart[i]));
    if(mv < min_map_value[i] || mv > max_map_value[i])  {
      status=-1;
      break;
    }
    mv = mv * weights[i];
    result += mv;
  }
  //return af::tiny<SiteFloatType, 2> (status, result);

//  MapFloatType status = 1;
//  for(std::size_t i_site=0;i_site<reference_selection.size();i_site++) {
//    std::size_t i = reference_selection[i_site];
//    MapFloatType mv = eight_point_interpolation(
//      density_map,
//      unit_cell.fractionalize(sites_cart[i]));
//    mv = mv * weights[i];
//    if(mv < reference_values[i]) status=-1;
//    //if(mv < reference_values[i] && mv > map_value_CB/3*weights[i]) status=-1;
//    //if(mv > reference_values[i]) {
//    //  reference_values[i] = mv;
//    //}
//  }


  //MapFloatType diff = 0;
  //for(std::size_t i=0; i<bonded_pairs.size(); i++) {
  //  af::tiny<std::size_t, 2> p = bonded_pairs[i];
  //  MapFloatType mv1 = eight_point_interpolation(
  //    density_map,
  //    unit_cell.fractionalize(sites_cart[p[0]]));
  //  MapFloatType mv2 = eight_point_interpolation(
  //    density_map,
  //    unit_cell.fractionalize(sites_cart[p[1]]));
  //  mv1 = std::abs(mv1);
  //  mv2 = std::abs(mv2);
  //  MapFloatType result;
  //  if(mv1>mv2) result = mv1/mv2;
  //  else        result = mv2/mv1;
  //  if(result > 5) {
  //    status = -1;
  //    std::cout<<result<<" "<<mv1<<" "<<mv2<<std::endl;
  //  }
  //  //mv1 = std::abs(mv1)*weights[p[0]];
  //  //mv2 = std::abs(mv2)*weights[p[1]];
  //  //diff += (std::abs(mv1-mv2));// *
  //  //        // std::abs(std::max(mv1*weights[p[0]], mv2*weights[p[1]])));
  //}
  return af::tiny<SiteFloatType, 2> (status, result);
}

// Keep for now!
//
//template <
//  typename MapFloatType,
//  typename SiteFloatType>
//MapFloatType
//target(
//  uctbx::unit_cell const& unit_cell,
//  af::const_ref<MapFloatType, af::c_grid_padded<3> > const& density_map,
//  af::const_ref<scitbx::vec3<SiteFloatType> > const& sites_cart,
//  af::const_ref<std::size_t> const& selection,
//  af::const_ref<MapFloatType > const& reference_values)
//{
//  CCTBX_ASSERT(sites_cart.size() == reference_values.size());
//  CCTBX_ASSERT(sites_cart.size() >= af::max(selection));
//  MapFloatType result = 0;
//  for(std::size_t i_site=0;i_site<sites_cart.size();i_site++) {
//    if(selection[i_site]) {
//      MapFloatType mv = eight_point_interpolation(
//        density_map,
//        unit_cell.fractionalize(sites_cart[i_site]));
//      result += mv;
//      //MapFloatType penalty = 0;
//      MapFloatType diff = mv - reference_values[i_site];
//      //if(diff<0) result += diff;
//    }
//  }
//  return result;
//}

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
target_with_adjacent_similarity(
  uctbx::unit_cell const& unit_cell,
  af::const_ref<MapFloatType, af::c_grid_padded<3> > const& density_map,
  af::const_ref<scitbx::vec3<SiteFloatType> > const& sites_cart,
  af::const_ref<std::size_t> const& selection,
  af::const_ref<SiteFloatType> const& weights)
{
  MapFloatType result = 0;
  af::shared<MapFloatType> vals;
  for(std::size_t i_site=0;i_site<selection.size();i_site++) {
    MapFloatType mv = eight_point_interpolation(
      density_map,
      unit_cell.fractionalize(sites_cart[selection[i_site]]));
    if(weights[i_site] != 0.) mv = mv/weights[i_site];
    result += mv;
    vals.push_back(mv);
  }
  // Penalty for adjacent values being too dissimilar
  for(std::size_t i=0;i<vals.size();i++) {
    if(i>=0 && i+1<vals.size()) {
      SiteFloatType v1 = vals[i];
      SiteFloatType v2 = vals[i+1];
      result -= std::abs(vals[i]-vals[i+1]);
    }
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
MapFloatType
target(
  af::const_ref<MapFloatType, af::c_grid_padded<3> > const& density_map,
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
