#ifndef CCTBX_MAPTBX_UTILS_H
#define CCTBX_MAPTBX_UTILS_H

#include <cstddef>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <cctbx/uctbx.h>
#include <scitbx/math/utils.h>
#include <scitbx/math/modulo.h>
#include <scitbx/array_family/sort.h>
#include <scitbx/random.h>
#include <scitbx/math/linear_correlation.h>

#include <cctbx/uctbx.h>

#include <cctbx/maptbx/interpolation.h> //indirect import?

#include <math.h>

namespace cctbx { namespace maptbx {

  //! Miller index element corresponding to 1-dimensional array index.
  template <typename IntegerType>
  inline
  IntegerType
  ih_as_h(IntegerType ih, std::size_t n_real)
  {
    if (ih <= n_real/2) return ih;
    return ih - n_real;
  }

  //! 1-dimensional array index corresponding to Miller index element.
  /*! Returns -1 if h is out of range (see code).
   */
  template <typename IntegerType>
  inline
  IntegerType
  h_as_ih_exact(IntegerType h, IntegerType n_complex, bool positive_only)
  {
    if (positive_only) {
      if (0 > h || h >= n_complex) return -1;
    }
    else {
      IntegerType m = (n_complex - 1) / 2;
      if (-m > h || h > m) return -1;
      else if (h < 0) return h + n_complex;
    }
    return h;
  }

  //! 3-dimensional array indices corresponding to Miller index.
  /*! Result is -1 for out-of-range elements (see code).
   */
  template <typename IndexTypeN>
  af::int3
  h_as_ih_exact_array(bool anomalous_flag,
                      miller::index<> const& h,
                      IndexTypeN const& n_complex)
  {
    af::int3 ih;
    bool positive_only[] = {false, false, !anomalous_flag};
    for(std::size_t i=0;i<3;i++) {
      ih[i] = h_as_ih_exact(h[i], n_complex[i], positive_only[i]);
    }
    return ih;
  }

  //! 1-dimensional array index corresponding to Miller index element.
  /*! Applies modulus operation (see code).
   */
  template <typename IntegerType>
  inline
  IntegerType
  h_as_ih_mod(IntegerType h, IntegerType const& n_real)
  {
    h %= n_real;
    if (h < 0) h += n_real;
    return h;
  }

  //! 3-dimensional array indices corresponding to Miller index.
  /*! Applies modulus operation (see code).
      <p>
      See also: structure_factors::to_map
   */
  template <typename IndexTypeN>
  inline
  IndexTypeN
  h_as_ih_mod_array(miller::index<> const& h, IndexTypeN const& n_real)
  {
    IndexTypeN ih;
    for(std::size_t i=0;i<3;i++) {
      ih[i] = h_as_ih_mod(h[i], n_real[i]);
    }
    return ih;
  }

template <typename DataType>
void hoppe_gassman_modification2(af::ref<DataType, af::c_grid<3> > map_data,
       int n_iterations)
/* A modified version of rho->3*rho^2-2*rho^3 modification.
   Acta Cryst. (1968). B24, 97-107
   Acta Cryst. (1975). A31, 388-389
   Acta Cryst. (1979). B35, 1776-1785
*/
{
  int nx = map_data.accessor()[0];
  int ny = map_data.accessor()[1];
  int nz = map_data.accessor()[2];
  for(int iter = 0; iter < n_iterations; iter++) {
    for(int i = 0; i < nx; i++) {
      for(int j = 0; j < ny; j++) {
        for(int k = 0; k < nz; k++) {
           DataType rho = map_data(i,j,k);
           if(rho<0) map_data(i,j,k) = 0;
           if(rho >=0 && rho<=1) {
             DataType rho_sq = rho*rho;
             map_data(i,j,k) = 3*rho_sq - 2*rho*rho_sq;
    }}}}
  }
}

template <typename DataType>
void hoppe_gassman_modification(af::ref<DataType, af::c_grid<3> > map_data,
       DataType mean_scale, int n_iterations)
/* A modified version of rho->3*rho^2-2*rho^3 modification.
   Acta Cryst. (1968). B24, 97-107
   Acta Cryst. (1975). A31, 388-389
   Acta Cryst. (1979). B35, 1776-1785
*/
{
  int nx = map_data.accessor()[0];
  int ny = map_data.accessor()[1];
  int nz = map_data.accessor()[2];
  for(int iter = 0; iter < n_iterations; iter++) {
    DataType rho_mean = 0;
    DataType rho_max = 0;
    int cntr = 0;
    for(int i = 0; i < nx; i++) {
      for(int j = 0; j < ny; j++) {
        for(int k = 0; k < nz; k++) {
          DataType rho = map_data(i,j,k);
          if(rho>0) {
            rho_mean += rho;
            cntr += 1;
            if(rho>rho_max) rho_max = rho;
    }}}}
    if(cntr != 0) rho_mean /= cntr;
    DataType rho_ms = rho_mean*mean_scale;
    if(rho_max!=0) {
      for(int i = 0; i < nx; i++) {
        for(int j = 0; j < ny; j++) {
          for(int k = 0; k < nz; k++) {
             DataType rho = map_data(i,j,k);
             if(rho > rho_ms) rho = rho_ms;
             CCTBX_ASSERT(rho<=rho_max);
             rho /= rho_max;
             if(rho<0) map_data(i,j,k) = 0;
             else {
               DataType rho_sq = rho*rho;
               map_data(i,j,k) = 3*rho_sq - 2*rho*rho_sq;
    }}}}}
  }
}

template <typename ComplexType, typename FloatType>
af::shared<FloatType> cc_complex_complex(
  af::const_ref<ComplexType> const& f_1,
  af::const_ref<ComplexType> const& f_2,
  af::const_ref<FloatType> const& d_spacings,
  af::const_ref<FloatType> const& ss,
  af::const_ref<FloatType> const& d_mins,
  FloatType const& b_iso)
{
  CCTBX_ASSERT(f_1.size()==f_2.size());
  CCTBX_ASSERT(f_1.size()==d_spacings.size());
  CCTBX_ASSERT(f_1.size()==ss.size());
  af::shared<FloatType> num_term(ss.size());
  af::shared<FloatType> d2_term(ss.size());
  af::shared<FloatType> f_1_i_sq(ss.size());
  FloatType d1 = 0.;
  for(int i = 0; i < ss.size(); i++) {
    FloatType f_2_scaled = std::abs(f_2[i] * std::exp(-b_iso*ss[i]));
    FloatType f_1_abs = std::abs(f_1[i]);
    FloatType f_1_arg = std::arg(f_1[i]);
    FloatType f_2_arg = std::arg(f_2[i]);
    num_term[i] = f_1_abs*f_2_scaled*std::cos(f_2_arg-f_1_arg);
    d2_term[i] = f_2_scaled*f_2_scaled;
    f_1_i_sq[i] = f_1_abs*f_1_abs;
    d1 += f_1_i_sq[i];
  }

  FloatType cc_best=-1.;
  FloatType d_min_best=-1.;
  af::shared<FloatType> result;
  int ss_size = ss.size();
  for(int j = 0; j < d_mins.size(); j++) {
    FloatType d_min = d_mins[j];
    FloatType num = 0.;
    FloatType den = 0.;
    FloatType d2 = 0.;
    FloatType cc = 0.;
    for(int i = 0; i < ss_size; i++) {
      if(d_spacings[i]>d_min) {
        num += num_term[i];
        d2 += d2_term[i];
      }
    }
    cc = num / (std::sqrt(d1*d2));
    if(cc>cc_best) {
      cc_best = cc;
      d_min_best = d_min;
    }
  }
  result.push_back(d_min_best);
  result.push_back(cc_best);
  return result;
}

template <typename ComplexType, typename FloatType>
FloatType cc_complex_complex(
  af::const_ref<ComplexType> const& f_1,
  af::const_ref<ComplexType> const& f_2)
{
  CCTBX_ASSERT(f_1.size()==f_2.size());
  af::shared<FloatType> num_term(f_1.size());
  af::shared<FloatType> d2_term(f_1.size());
  af::shared<FloatType> f_1_i_sq(f_1.size());
  for(int i = 0; i < f_1.size(); i++) {
    FloatType f_2_abs = std::abs(f_2[i]);
    FloatType f_1_abs = std::abs(f_1[i]);
    FloatType f_1_arg = std::arg(f_1[i]);
    FloatType f_2_arg = std::arg(f_2[i]);
    num_term[i] = f_1_abs*f_2_abs*std::cos(f_2_arg-f_1_arg);
    d2_term[i] = f_2_abs*f_2_abs;
    f_1_i_sq[i] = f_1_abs*f_1_abs;
  }
  FloatType num = 0.;
  FloatType d1 = 0.;
  FloatType d2 = 0.;
  for(int i = 0; i < f_1.size(); i++) {
    num += num_term[i];
    d2 += d2_term[i];
    d1 += f_1_i_sq[i];
  }
  FloatType cc = num / (std::sqrt(d1*d2));
  return cc;
}

/*
Acta Cryst. (2014). D70, 2593-2606
Metrics for comparison of crystallographic maps
A. Urzhumtsev, P. V. Afonine, V. Y. Lunin, T. C. Terwilliger and P. D. Adams
*/
template <typename DataType>
DataType cc_peak(
  af::const_ref<DataType, af::c_grid<3> > const& map_1,
  af::const_ref<DataType, af::c_grid<3> > const& map_2,
  DataType const& cutoff)
{
  af::c_grid<3> a1 = map_1.accessor();
  af::c_grid<3> a2 = map_2.accessor();
  for(int i = 0; i < 3; i++) CCTBX_ASSERT(a1[i]==a2[i]);
  af::shared<DataType> m1;
  af::shared<DataType> m2;
  for(int i = 0; i < a1[0]; i++) {
    for(int j = 0; j < a1[1]; j++) {
      for(int k = 0; k < a1[2]; k++) {
        DataType m1_ = map_1(i,j,k);
        DataType m2_ = map_2(i,j,k);
        if(m1_>=cutoff && m2_>=cutoff) {
          m1.push_back(m1_);
          m2.push_back(m2_);
        }
        else if(m1_<cutoff && m2_>=cutoff) {
          m1.push_back(cutoff);
          m2.push_back(m2_);
        }
        else if(m1_>=cutoff && m2_<cutoff) {
          m1.push_back(m1_);
          m2.push_back(cutoff);
        }
  }}}
  return scitbx::math::linear_correlation<>(
    m1.ref(), m2.ref()).coefficient();;
}

template <typename DataType>
DataType map_sum_at_sites_frac(
  af::const_ref<DataType, af::c_grid<3> > const& map_data,
  af::const_ref<scitbx::vec3<DataType> > const& sites_frac)
{
  DataType result = 0;
  for(int i = 0; i < sites_frac.size(); i++) {
    result += tricubic_interpolation(map_data, sites_frac[i]);
  }
  return result;
}

template <typename DataType>
af::versa<DataType, af::c_grid<3> > set_box_copy_inside(
  DataType const& value,
  af::ref<DataType, af::c_grid<3> > map_data_to,
  af::tiny<int, 3> const& start,
  af::tiny<int, 3> const& end)
{
  af::c_grid<3> a = map_data_to.accessor();
  for(int i = 0; i < 3; i++) {
    CCTBX_ASSERT(start[i]>=0 && start[i]<=a[i]);
    CCTBX_ASSERT(end[i]>=0   && end[i]<=a[i]);
  }
  // like set_box_copy except copies inside and sets value outside
  af::versa<DataType, af::c_grid<3> > result_map(a,
    af::init_functor_null<DataType>());
  af::ref<DataType, af::c_grid<3> > result_map_ref = result_map.ref();
  for(int i = 0; i < a[0]; i++) {
    for(int j = 0; j < a[1]; j++) {
      for(int k = 0; k < a[2]; k++) {
        if(i>=start[0]&&i<end[0] &&
           j>=start[1]&&j<end[1] &&
           k>=start[2]&&k<end[2]) {
          result_map_ref(i,j,k) = map_data_to(i,j,k);
        }
        else {
          result_map_ref(i,j,k) = value;
        }
  }}}
  return result_map;
}

template <typename DataType>
af::versa<DataType, af::c_grid<3> > set_box_copy(
  DataType const& value,
  af::ref<DataType, af::c_grid<3> > map_data_to,
  af::tiny<int, 3> const& start,
  af::tiny<int, 3> const& end)
{
  af::c_grid<3> a = map_data_to.accessor();
  for(int i = 0; i < 3; i++) {
    CCTBX_ASSERT(start[i]>=0 && start[i]<=a[i]);
    CCTBX_ASSERT(end[i]>=0   && end[i]<=a[i]);
  }
  af::versa<DataType, af::c_grid<3> > result_map(a,
    af::init_functor_null<DataType>());
  af::ref<DataType, af::c_grid<3> > result_map_ref = result_map.ref();
  for(int i = 0; i < a[0]; i++) {
    for(int j = 0; j < a[1]; j++) {
      for(int k = 0; k < a[2]; k++) {
        if(i>=start[0]&&i<end[0] &&
           j>=start[1]&&j<end[1] &&
           k>=start[2]&&k<end[2]) {
          result_map_ref(i,j,k) = value;
        }
        else {
          result_map_ref(i,j,k) = map_data_to(i,j,k);
        }
  }}}
  return result_map;
}

template <typename ComplexType, typename FloatType>
af::shared<ComplexType> fem_averaging_loop(
  af::const_ref<ComplexType> const& map_coefficients,
  af::const_ref<FloatType> const& r_factors,
  af::const_ref<FloatType> const& sigma_over_f_obs,
  FloatType const& random_scale,
  int const& random_seed,
  int const& n_cycles)
{
  CCTBX_ASSERT(n_cycles>0);
  CCTBX_ASSERT(r_factors.size()==sigma_over_f_obs.size());
  CCTBX_ASSERT(r_factors.size()==map_coefficients.size());
  af::shared<ComplexType> result(r_factors.size());
  for(int i = 0; i < result.size(); i++) result[i] = ComplexType(0,0);
  scitbx::random::mersenne_twister mt(random_seed);
  for(int j = 0; j < n_cycles; j++) {
    for(int i = 0; i < map_coefficients.size(); i++) {
      FloatType s1 = mt.random_double()*random_scale;
      FloatType s2 = mt.random_double()*random_scale;
      FloatType one_over_w = 1. + r_factors[i]*s1 + sigma_over_f_obs[i]*s2;
      CCTBX_ASSERT(one_over_w != 0);
      result[i] = result[i] + map_coefficients[i] / one_over_w;
  }}
  for(int i = 0; i < result.size(); i++) result[i] = result[i] * (1./n_cycles);
  return result;
}

template <typename FloatType>
cctbx::cartesian<>
center_of_mass(
  af::const_ref<FloatType, af::c_grid<3> > const& map_data,
  cctbx::uctbx::unit_cell const& unit_cell,
  FloatType const& cutoff)
{
  FloatType result = 0.;
  FloatType mass_sum = 0.;
  cctbx::cartesian<> num = cctbx::cartesian<> (0,0,0);
  af::c_grid<3> a = map_data.accessor();
  for (int i = 0; i < a[0]; i++) {
    for (int j = 0; j < a[1]; j++) {
      for (int k = 0; k < a[2]; k++) {
        FloatType m = map_data(i,j,k);
        if(m > cutoff) {
          cctbx::fractional<> grid_frac = cctbx::fractional<>(
            i/static_cast<FloatType>(a[0]),
            j/static_cast<FloatType>(a[1]),
            k/static_cast<FloatType>(a[2]));
          num += (unit_cell.orthogonalize(grid_frac)*m);
          mass_sum += m;
  }}}}
  CCTBX_ASSERT(mass_sum != 0);
  return num/mass_sum;
}

template <typename DataType>
af::versa<DataType, af::c_grid<3> > conditional_solvent_region_filter(
  af::const_ref<DataType, af::c_grid_padded<3> > const& bulk_solvent_mask,
  af::const_ref<DataType, af::c_grid<3> > const& map_data,
  DataType const& threshold)
{
  af::tiny<int, 3> a2 = bulk_solvent_mask.accessor().all();
  af::c_grid<3> a1 = map_data.accessor();
  for(int i = 0; i < 3; i++) CCTBX_ASSERT(a1[i]==a2[i]);
  af::versa<DataType, af::c_grid<3> > result_map(a1,
    af::init_functor_null<DataType>());
  af::ref<DataType, af::c_grid<3> > result_map_ref = result_map.ref();
  for(int i = 0; i < a1[0]; i++) {
    for(int j = 0; j < a1[1]; j++) {
      for(int k = 0; k < a1[2]; k++) {
        DataType map_  = map_data(i,j,k);
        if(bulk_solvent_mask(i,j,k)==0) {
          result_map_ref(i,j,k) = 1;
        }
        else {
          if(map_<threshold) {
            result_map_ref(i,j,k) = 0;
          }
          else {
            result_map_ref(i,j,k) = 1;
          }
        }
  }}}
  return result_map;
}

template <typename DataType1, typename DataType2>
af::versa<DataType2, af::c_grid<3> > update_f_part1_helper(
  af::const_ref<DataType1, af::c_grid_padded<3> > const& connectivity_map,
  af::const_ref<DataType2, af::c_grid<3> > const& map_data,
  int const& region_id)
{
  af::tiny<int, 3> a2 = connectivity_map.accessor().all();
  af::c_grid<3> a1 = map_data.accessor();
  for(int i = 0; i < 3; i++) CCTBX_ASSERT(a1[i]==a2[i]);
  af::versa<DataType2, af::c_grid<3> > result_map(a1,
    af::init_functor_null<DataType2>());
  af::ref<DataType2, af::c_grid<3> > result_map_ref = result_map.ref();
  for(int i = 0; i < a1[0]; i++) {
    for(int j = 0; j < a1[1]; j++) {
      for(int k = 0; k < a1[2]; k++) {
        DataType1 cm = connectivity_map(i,j,k);
        DataType1 md = map_data(i,j,k);
        if(cm==region_id) {
          result_map_ref(i,j,k) = -1.*md;
        }
        else {
          result_map_ref(i,j,k) = 0;
        }
  }}}
  return result_map;
}

template <typename DataType>
void set_box(
  DataType const& value,
  af::ref<DataType, af::c_grid<3> > map_data_to,
  af::tiny<int, 3> const& start,
  af::tiny<int, 3> const& end)
{
  af::c_grid<3> a = map_data_to.accessor();
  af::tiny<int, 3> startmod, endmod;
  af::shared<DataType> gridstartx, gridstarty, gridstartz;
  af::shared<DataType> gridendx, gridendy, gridendz;
  // make sure that grid points are not further apart than unit cell length
  // and that start is strictly smaller than end
  for(int i = 0; i < 3; i++) {
    CCTBX_ASSERT((end[i] - start[i]) <= a[i]);
    CCTBX_ASSERT(end[i] > start[i]);
  }
  for(int i = 0; i < 3; i++) {
    startmod[i] = scitbx::math::mod_positive(start[i], static_cast<int>(a[i]));
    endmod[i] = scitbx::math::mod_positive(end[i], static_cast<int>(a[i]));
    if (endmod[i] == 0) endmod[i] = a[i];
  }
  //
  gridstartx.push_back(startmod[0]);
  gridendx.push_back(endmod[0]);
  if (startmod[0] > endmod[0]) {
    gridstartx.insert(gridstartx.begin(), 0);
    gridendx.push_back(a[0]);
  }
  gridstarty.push_back(startmod[1]);
  gridendy.push_back(endmod[1]);
  if (startmod[1] > endmod[1]) {
    gridstarty.insert(gridstarty.begin(), 0);
    gridendy.push_back(a[1]);
  }
  gridstartz.push_back(startmod[2]);
  gridendz.push_back(endmod[2]);
  if (startmod[2] > endmod[2]) {
    gridstartz.insert(gridstartz.begin(), 0);
    gridendz.push_back(a[2]);
  }
  for (int l = 0; l < gridstartx.size(); l++) {
    for (int m = 0; m < gridstarty.size(); m++) {
      for (int n = 0; n < gridstartz.size(); n++) {
        for (int i = gridstartx[l]; i < gridendx[l]; i++) {
          for (int j = gridstarty[m]; j < gridendy[m]; j++) {
            for (int k = gridstartz[n]; k < gridendz[n]; k++) {
              map_data_to(i,j,k) = value;
            }
          }
        }
      }
    }
  }
}

template <typename DataType>
void set_box(
  af::const_ref<DataType, af::c_grid<3> > const& map_data_from,
  af::ref<DataType, af::c_grid<3> > map_data_to,
  af::tiny<int, 3> const& start,
  af::tiny<int, 3> const& end)
{
  af::c_grid<3> a = map_data_to.accessor();
  //
  // PVA: I need to find out why this occasionally crashes here:
  //
  //for(int i = 0; i < 3; i++) {
  //  CCTBX_ASSERT((end[i] - start[i]) <= a[i]);
  //  CCTBX_ASSERT(end[i] > start[i]);
  //}
  for (int i = start[0]; i < end[0]; i++) {
    int ii = scitbx::math::mod_positive(i, static_cast<int>(a[0]));
    for (int j = start[1]; j < end[1]; j++) {
      int jj = scitbx::math::mod_positive(j, static_cast<int>(a[1]));
      for (int k = start[2]; k < end[2]; k++) {
        int kk = scitbx::math::mod_positive(k, static_cast<int>(a[2]));
        int p = i-start[0];
        int q = j-start[1];
        int r = k-start[2];
        map_data_to(ii,jj,kk) = map_data_from(p,q,r);
      }
    }
  }
}

template <typename DataType>
void copy_box(
  af::const_ref<DataType, af::c_grid<3> > const& map_data_from,
  af::ref<DataType, af::c_grid<3> > map_data_to,
  af::tiny<int, 3> const& start,
  af::tiny<int, 3> const& end)
{
  af::c_grid<3> a1 = map_data_to.accessor();
  af::c_grid<3> a2 = map_data_from.accessor();
  for(int i = 0; i < 3; i++) {
    CCTBX_ASSERT(a1[i]==a2[i]);
    CCTBX_ASSERT(start[i]>=0 && start[i]<=a1[i]);
    CCTBX_ASSERT(end[i]>=0   && end[i]<=a1[i]);
  }
  for (int i = start[0]; i < end[0]; i++) {
    for (int j = start[1]; j < end[1]; j++) {
      for (int k = start[2]; k < end[2]; k++) {
        map_data_to(i,j,k) = map_data_from(i,j,k);
      }
    }
  }
}

template <typename DataType>
DataType
_one_map_box_average(
  af::ref<DataType, af::c_grid<3> > map_data,
  int const& index_span,
  int const& lx,
  int const& ly,
  int const& lz)
{
  af::c_grid<3> a = map_data.accessor();
  DataType rho = 0.0;
  int counter = 0;
  for (int i = lx-index_span; i <= lx+index_span; i++) {
    int mx = i;
    if(i<0 || i>=a[0]) {
      mx = scitbx::math::mod_positive(i, static_cast<int>(a[0]));
    }
    for (int j = ly-index_span; j <= ly+index_span; j++) {
      int my = j;
      if(j<0 || j>=a[1]) {
        my = scitbx::math::mod_positive(j, static_cast<int>(a[1]));
      }
      for (int k = lz-index_span; k <= lz+index_span; k++) {
        int mz = k;
        if(k<0 || k>=a[2]) {
          mz = scitbx::math::mod_positive(k, static_cast<int>(a[2]));
        }
        rho += map_data(mx,my,mz);
        counter += 1;
  }}}
  return rho / counter;
}

template <typename DataType>
void
map_box_average(
  af::ref<DataType, af::c_grid<3> > map_data,
  int const& index_span)
{
  int nx = map_data.accessor()[0];
  int ny = map_data.accessor()[1];
  int nz = map_data.accessor()[2];
  for (int lx = 0; lx < nx; lx++) {
    for (int ly = 0; ly < ny; ly++) {
      for (int lz = 0; lz < nz; lz++) {
        map_data(lx,ly,lz) = _one_map_box_average(
          map_data, index_span, lx, ly, lz);
  }}}
}

template <typename DataType>
void
map_box_average(
  af::ref<DataType, af::c_grid<3> > map_data,
  int const& index_span,
  DataType const& threshold)
{
  int nx = map_data.accessor()[0];
  int ny = map_data.accessor()[1];
  int nz = map_data.accessor()[2];
  for (int lx = 0; lx < nx; lx++) {
    for (int ly = 0; ly < ny; ly++) {
      for (int lz = 0; lz < nz; lz++) {
        if(std::abs(map_data(lx,ly,lz))<threshold) {
          map_data(lx,ly,lz) = _one_map_box_average(
            map_data, index_span, lx, ly, lz);
        }
  }}}
}

template <typename DataType>
int nint(DataType x)
{
  // Analogue of FORTRAN NINT(X)
  return int(std::ceil(x+0.5)-(std::fmod(x*0.5+0.25,1.0)!=0));
}

template <typename DataType>
void set_box_with_symmetry(
  af::const_ref<DataType, af::c_grid<3> > const& map_data_from,
  af::ref<DataType, af::c_grid<3> > map_data_to,
  af::tiny<int, 3> const& start,
  af::tiny<int, 3> const& end,
  cctbx::uctbx::unit_cell const& unit_cell,
  af::shared<scitbx::mat3<double> > const& rotation_matrices,
  af::shared<scitbx::vec3<double> > const& translation_vectors)
{
  /*
  This is a specialized function. It takes two maps: a larger empty map
  (map_data_to) and smaller map (map_data_from).
  Parameter start are the coordinates of the origin of map_data_from wrt origin
  of map_data_to.
  Parameters start, end, unit_cell are related to map_data_to.
  */
  af::c_grid<3> a = map_data_to.accessor();
  af::c_grid<3> b = map_data_from.accessor();
  for (int i = start[0]; i <= end[0]; i++) {
    for (int j = start[1]; j <= end[1]; j++) {
      for (int k = start[2]; k <= end[2]; k++) {
        // position in map_data_from
        int p = scitbx::math::mod_positive(i-start[0], static_cast<int>(b[0]));
        int q = scitbx::math::mod_positive(j-start[1], static_cast<int>(b[1]));
        int r = scitbx::math::mod_positive(k-start[2], static_cast<int>(b[2]));
        // fractional coordinates of i,j,k
        cctbx::fractional<> grid_node_frac = cctbx::fractional<>(
          i/static_cast<double>(a[0]),
          j/static_cast<double>(a[1]),
          k/static_cast<double>(a[2]));
        for (int o = 0; o < rotation_matrices.size(); o++) {
          scitbx::mat3<double> rm = rotation_matrices[o];
          scitbx::vec3<double> tv = translation_vectors[o];
          cctbx::fractional<> grid_node_frac_rt = rm * grid_node_frac + tv;
          // position of rotated+translated point in original map
          int ii = scitbx::math::mod_positive(
            nint(grid_node_frac_rt[0]*a[0]), static_cast<int>(a[0]));
          int jj = scitbx::math::mod_positive(
            nint(grid_node_frac_rt[1]*a[1]), static_cast<int>(a[1]));
          int kk = scitbx::math::mod_positive(
            nint(grid_node_frac_rt[2]*a[2]), static_cast<int>(a[2]));
          // Using max avoids overwriting set non-zero values with zeros from the box
         double mv_from = map_data_from(p,q,r);
         double mv_to   = map_data_to(ii,jj,kk);
         if(std::abs(mv_to)<1.e-6) {
           map_data_to(ii,jj,kk) = mv_from;
         }
         else {
           if(std::abs(mv_from)>1.e-6) map_data_to(ii,jj,kk) = (mv_from+mv_to)/2;
         }
        }
      }
    }
  }
  map_box_average(map_data_to, 1, 1.e-6);
}

template <typename DataType1, typename DataType2>
void combine_1(
       af::ref<DataType2, af::c_grid<3> > map_data,
       af::ref<DataType1, af::c_grid<3> > diff_map)
{
  af::tiny<int, 3> a = map_data.accessor();
  for(int i = 0; i < a[0]; i++) {
    for(int j = 0; j < a[1]; j++) {
      for(int k = 0; k < a[2]; k++) {
         DataType1 m1 = map_data(i,j,k);
         DataType2 m2 = diff_map(i,j,k);
         if(m1<=1.0) map_data(i,j,k) = m1+m2;
  }}}
}

template <typename DataType1, typename DataType2>
void truncate_special(
       af::ref<DataType2, af::c_grid<3> > mask,
       af::ref<DataType1, af::c_grid<3> > map_data)
{
  af::tiny<int, 3> a = map_data.accessor();
  for(int i = 0; i < a[0]; i++) {
    for(int j = 0; j < a[1]; j++) {
      for(int k = 0; k < a[2]; k++) {
         DataType1 m1 = map_data(i,j,k);
         DataType2 m2 = mask(i,j,k);
         if(m2<1 && m1<=0.5) map_data(i,j,k) = 0;
  }}}
}

template <typename DataType>
void truncate_between_min_max(
       af::ref<DataType, af::c_grid<3> > map_data,
       DataType const& min,
       DataType const& max)
{
  af::tiny<int, 3> a = map_data.accessor();
  for(int i = 0; i < a[0]; i++) {
    for(int j = 0; j < a[1]; j++) {
      for(int k = 0; k < a[2]; k++) {
         DataType md = map_data(i,j,k);
         if(md > min && md < max) map_data(i,j,k) = 0;
  }}}
}

template <typename DataType>
void binarize(
       af::ref<DataType, af::c_grid<3> > map_data,
       DataType const& threshold,
       DataType const& substitute_value_below,
       DataType const& substitute_value_above)
{
  af::tiny<int, 3> a = map_data.accessor();
  for(int i = 0; i < a[0]; i++) {
    for(int j = 0; j < a[1]; j++) {
      for(int k = 0; k < a[2]; k++) {
         DataType* md = &map_data(i,j,k);
         if(*md<threshold) *md = substitute_value_below;
         else              *md = substitute_value_above;
  }}}
}

template <typename DataType>
void truncate(
       af::ref<DataType, af::c_grid<3> > map_data,
       DataType const& standard_deviation,
       DataType const& by_sigma_less_than,
       DataType const& scale_by,
       DataType const& set_value)
{
  af::tiny<int, 3> a = map_data.accessor();
  for(int i = 0; i < a[0]; i++) {
    for(int j = 0; j < a[1]; j++) {
      for(int k = 0; k < a[2]; k++) {
         DataType md = map_data(i,j,k);
         if(md/standard_deviation > by_sigma_less_than) md = set_value*scale_by;
         else md = md * scale_by;
         map_data(i,j,k)=md;
  }}}
}

template <typename DataType>
void intersection(
       af::ref<DataType, af::c_grid<3> > map_data_1,
       af::ref<DataType, af::c_grid<3> > map_data_2,
       af::ref<DataType> const& thresholds,
       bool average)
{
  af::tiny<int, 3> a1 = map_data_1.accessor();
  af::tiny<int, 3> a2 = map_data_2.accessor();
  for(int i = 0; i < 3; i++) CCTBX_ASSERT(a1[i]==a2[i]);
  for(int i = 0; i < a1[0]; i++) {
    for(int j = 0; j < a1[1]; j++) {
      for(int k = 0; k < a1[2]; k++) {
        double rho1 = map_data_1(i,j,k);
        double rho2 = map_data_2(i,j,k);
        for(int m = 0; m < thresholds.size(); m++) {
          DataType threshold = thresholds[m];
          bool c1 = rho1>threshold && rho2<threshold;
          bool c2 = rho2>threshold && rho1<threshold;
          if(c1 || c2) {
            map_data_1(i,j,k)=0;
            map_data_2(i,j,k)=0;
        }
        if(average) {
          DataType rho_ave = (map_data_1(i,j,k)+map_data_2(i,j,k))/2.;
          map_data_1(i,j,k) = rho_ave;
          map_data_2(i,j,k) = rho_ave;
        }
  }}}}
}

template <typename DataType>
void convert_to_non_negative(
       af::ref<DataType, af::c_grid<3> > map_data,
       DataType substitute_value)
{
  int nx = map_data.accessor()[0];
  int ny = map_data.accessor()[1];
  int nz = map_data.accessor()[2];
  double rho_max = af::max(map_data);
  for(int i = 0; i < nx; i++) {
    for(int j = 0; j < ny; j++) {
      for(int k = 0; k < nz; k++) {
         double rho = map_data(i,j,k);
         if(rho<0) map_data(i,j,k) = substitute_value;
  }}}
}

template <typename DataType>
void flexible_boundary_mask(
       af::ref<DataType, af::c_grid<3> > map_data,
       af::ref<DataType, af::c_grid<3> > mask_data)
{
  int nx = map_data.accessor()[0];
  int ny = map_data.accessor()[1];
  int nz = map_data.accessor()[2];
  for(int i = 0; i < nx; i++) {
    for(int j = 0; j < ny; j++) {
      for(int k = 0; k < nz; k++) {
         double r = map_data(i,j,k);
         double m = mask_data(i,j,k);
         mask_data(i,j,k)=std::max(m-r, 0.0);
  }}}
}

template <typename DataType>
void reset(
       af::ref<DataType, af::c_grid<3> > map_data,
       DataType substitute_value,
       DataType less_than_threshold,
       DataType greater_than_threshold,
       bool use_and)
{
  int nx = map_data.accessor()[0];
  int ny = map_data.accessor()[1];
  int nz = map_data.accessor()[2];
  for(int i = 0; i < nx; i++) {
    for(int j = 0; j < ny; j++) {
      for(int k = 0; k < nz; k++) {
         double rho = map_data(i,j,k);
         if(use_and) {
           if(rho<less_than_threshold && rho>greater_than_threshold) {
             map_data(i,j,k) = substitute_value;
           }
         }
         else {
           if(rho<less_than_threshold || rho>greater_than_threshold) {
             map_data(i,j,k) = substitute_value;
           }
         }
  }}}
}

template <typename DataType>
void resample(
       af::const_ref<DataType, af::c_grid<3> > const& map_data,
       af::ref<DataType, af::c_grid<3> > map_data_new,
       cctbx::uctbx::unit_cell const& unit_cell)
{
  af::tiny<int, 3> const& n_real = map_data_new.accessor();
  af::tiny<DataType, 6> ucp = unit_cell.parameters();
  DataType sx = ucp[0]/n_real[0];
  DataType sy = ucp[1]/n_real[1];
  DataType sz = ucp[2]/n_real[2];
  for(int i = 0; i < n_real[0]; i++) {
    DataType x = sx*i;
    for(int j = 0; j < n_real[1]; j++) {
      DataType y = sy*j;
      for(int k = 0; k < n_real[2]; k++) {
        cctbx::cartesian<>  site_cart = cctbx::cartesian<>(x, y, sz*k);
        cctbx::fractional<> site_frac = unit_cell.fractionalize(site_cart);
        map_data_new(i,j,k) = tricubic_interpolation(map_data, site_frac);
  }}}
}

template <typename DataType>
void negate_selected_in_place(
       af::ref<DataType, af::c_grid<3> > map_data,
       af::const_ref<bool, af::c_grid<3> > const& selection)
{
  int nx = map_data.accessor()[0];
  int ny = map_data.accessor()[1];
  int nz = map_data.accessor()[2];
  for(int i = 0; i < nx; i++) {
    for(int j = 0; j < ny; j++) {
      for(int k = 0; k < nz; k++) {
        if(selection(i,j,k)) {
          map_data(i,j,k) = -1.*map_data(i,j,k);
        }
  }}}
}

template <typename DataType, typename GridType>
af::versa<DataType, GridType>
negate_selected_in_place(
  af::const_ref<DataType, GridType> const& map_data,
  std::vector<unsigned> const& selection)
{
  af::c_grid_padded<3> a = map_data.accessor();
  int nx = map_data.accessor().all()[0];
  int ny = map_data.accessor().all()[1];
  int nz = map_data.accessor().all()[2];
  af::versa<DataType, af::c_grid_padded<3> > result_map(a,
    af::init_functor_null<DataType>());
  af::ref<DataType, af::c_grid_padded<3> > result_map_ref = result_map.ref();
  for(int i = 0; i < nx; i++) {
    for(int j = 0; j < ny; j++) {
      for(int k = 0; k < nz; k++) {
          result_map_ref(i,j,k) = map_data(i,j,k);
    }}}
  for(int i = 0; i < selection.size(); i++) {
     result_map[selection[i]] = -1.*map_data[selection[i]];
  }
  return result_map;
}

template <typename DataType>
void
median_filter(
  af::ref<DataType, af::c_grid<3> > map_data,
  int const& index_span)
{
  af::shared<DataType> box;
  int nx = map_data.accessor()[0];
  int ny = map_data.accessor()[1];
  int nz = map_data.accessor()[2];
  for (int lx = 0; lx < nx; lx++) {
    for (int ly = 0; ly < ny; ly++) {
      for (int lz = 0; lz < nz; lz++) {
        box.resize(0,0);
        for (int i = lx-index_span; i <= lx+index_span; i++) {
          for (int j = ly-index_span; j <= ly+index_span; j++) {
            for (int k = lz-index_span; k <= lz+index_span; k++) {
              int mx = scitbx::math::mod_positive(i, nx);
              int my = scitbx::math::mod_positive(j, ny);
              int mz = scitbx::math::mod_positive(k, nz);
              DataType rho = map_data(mx,my,mz);
              box.push_back(rho);
        }}}
        af::shared<std::size_t> permut;
        permut = af::sort_permutation(box.const_ref(),true);
        DataType median = box[permut[int(permut.size()/2)]];
        map_data(lx,ly,lz) = median;
  }}}
}

template <typename DataType>
void
kuwahara_filter(
  af::ref<DataType, af::c_grid<3> > map_data,
  int const& index_span)
{
  af::tiny<af::shared<DataType>, 8> octants;
  af::tiny<DataType, 8> mean;
  af::tiny<DataType, 8> variance;
  int nx = map_data.accessor()[0];
  int ny = map_data.accessor()[1];
  int nz = map_data.accessor()[2];
  for(int lx = 0; lx < nx; lx++) {
    for(int ly = 0; ly < ny; ly++) {
      for(int lz = 0; lz < nz; lz++) {
        // extract octants
        for(int o = 0; o < 8; o++) {
          af::shared<DataType> tmp;
          octants[o] = tmp;
        }
        for(int i = -index_span; i <= index_span; i++) {
          for(int j = -index_span; j <= index_span; j++) {
            for(int k = -index_span; k <= index_span; k++) {
              int mx = scitbx::math::mod_positive(lx+i, nx);
              int my = scitbx::math::mod_positive(ly+j, ny);
              int mz = scitbx::math::mod_positive(lz+k, nz);
              DataType rho = map_data(mx,my,mz);
              if(i>=0&&j>=0&&k>=0) octants[0].push_back(rho); // octant #1
              if(i<=0&&j>=0&&k>=0) octants[1].push_back(rho); // octant #2
              if(i<=0&&j<=0&&k>=0) octants[2].push_back(rho); // octant #3
              if(i>=0&&j<=0&&k>=0) octants[3].push_back(rho); // octant #4
              if(i>=0&&j>=0&&k<=0) octants[4].push_back(rho); // octant #5
              if(i<=0&&j>=0&&k<=0) octants[5].push_back(rho); // octant #6
              if(i<=0&&j<=0&&k<=0) octants[6].push_back(rho); // octant #7
              if(i>=0&&j<=0&&k<=0) octants[7].push_back(rho); // octant #8
        }}}
        // compute mean and variance in each octant
        mean.fill(0);
        variance.fill(0);
        for(int o = 0; o < 8; o++) {
          af::shared<DataType> oct = octants[o];
          CCTBX_ASSERT(oct.size()!=0);
          // mean
          DataType m = 0;
          for(int p = 0; p < oct.size(); p++) m += oct[p];
          mean[o] = m/oct.size();
          // variance
          DataType v = 0;
          for(int p = 0; p < oct.size(); p++) {
            DataType delta = oct[p] - mean[o];
            v += (delta*delta);
          }
          variance[o] = v/oct.size();
        }
        // find mean corresponding to the lowest variance
        DataType mean_result = 0;
        DataType v_smallest = 1.e99;
        for(int o = 0; o < 8; o++) {
          if(variance[o] < v_smallest) {
            v_smallest = variance[o];
            mean_result = mean[o];
          }
        }
        // set to central pixel
        map_data(lx,ly,lz) = mean_result;
  }}}
}

template <typename DataType>
void
remove_single_node_peaks(
  af::ref<DataType, af::c_grid<3> > map_data,
  af::ref<DataType, af::c_grid_padded<3> > mask_data,
  DataType const& cutoff,
  int const& index_span)
{
  int nx = map_data.accessor()[0];
  int ny = map_data.accessor()[1];
  int nz = map_data.accessor()[2];
  for (int lx = 0; lx < nx; lx++) {
    for (int ly = 0; ly < ny; ly++) {
      for (int lz = 0; lz < nz; lz++) {
          if(mask_data(lx,ly,lz)==0) continue;
          int counter = 0;
          for (int i = lx-index_span; i <= lx+index_span; i=i+2) {
            for (int j = ly-index_span; j <= ly+index_span; j=j+2) {
              for (int k = lz-index_span; k <= lz+index_span; k=k+2) {
                if(i!=lx || j!=ly || k!=lz) {
                  int mx = scitbx::math::mod_positive(i, nx);
                  int my = scitbx::math::mod_positive(j, ny);
                  int mz = scitbx::math::mod_positive(k, nz);
                  if(map_data(mx,my,mz)< cutoff) counter += 1;
                }
          }}}
          if(counter == 26) {
            map_data(lx,ly,lz) = 0;
          }
  }}}
}

template <typename DataType>
void
map_box_average(
  af::ref<DataType, af::c_grid<3> > map_data,
  cctbx::uctbx::unit_cell const& unit_cell,
  double const& radius)
{
  int nx = map_data.accessor()[0];
  int ny = map_data.accessor()[1];
  int nz = map_data.accessor()[2];
  DataType xrad = radius*unit_cell.reciprocal_parameters()[0]*nx;
  DataType yrad = radius*unit_cell.reciprocal_parameters()[1]*ny;
  DataType zrad = radius*unit_cell.reciprocal_parameters()[2]*nz;
  for (int lx = 0; lx < nx; lx++) {
    for (int ly = 0; ly < ny; ly++) {
      for (int lz = 0; lz < nz; lz++) {
        DataType rho = 0.0;
        int counter = 0;
        int x1box=scitbx::math::nearest_integer(static_cast<DataType>(lx)-xrad);
        int x2box=scitbx::math::nearest_integer(static_cast<DataType>(lx)+xrad);
        int y1box=scitbx::math::nearest_integer(static_cast<DataType>(ly)-yrad);
        int y2box=scitbx::math::nearest_integer(static_cast<DataType>(ly)+yrad);
        int z1box=scitbx::math::nearest_integer(static_cast<DataType>(lz)-zrad);
        int z2box=scitbx::math::nearest_integer(static_cast<DataType>(lz)+zrad);
        for (int kx = x1box; kx <= x2box; kx++) {
          for (int ky = y1box; ky <= y2box; ky++) {
            for (int kz = z1box; kz <= z2box; kz++) {
              int mx = scitbx::math::mod_positive(kx, nx);
              int my = scitbx::math::mod_positive(ky, ny);
              int mz = scitbx::math::mod_positive(kz, nz);
              rho += map_data(mx,my,mz);
              counter += 1;
        }}}
        map_data(lx,ly,lz) = rho / counter;
  }}}
}

class sphericity2 {
public:
  af::tiny<double, 3> rho_min_max_mean_;
  af::tiny<double, 3> ccs_min_max_mean_;
  sphericity2(
    af::const_ref<double, af::c_grid<3> > const& map_data,
    cctbx::cartesian<double> const& center_cart,
    af::const_ref<scitbx::vec3<double> > const& points_on_sphere_cart,
    cctbx::uctbx::unit_cell const& unit_cell)
  :
  rho_min_max_mean_(af::tiny<double, 3>(0,0,0)),
  ccs_min_max_mean_(af::tiny<double, 3>(0,0,0))
  {
    af::shared<af::shared<double> > vecs(points_on_sphere_cart.size());
    double rho_min=1.e+9, rho_max=-1.e+9, rho_mean=0.;
    for(int i=0; i<points_on_sphere_cart.size(); i++) {
      cctbx::cartesian<double>  posc = points_on_sphere_cart[i];
      cctbx::fractional<double> posf = unit_cell.fractionalize(posc);
      double rho = tricubic_interpolation(map_data, posf);
      if(rho<rho_min) rho_min = rho;
      if(rho>rho_max) rho_max = rho;
      rho_mean += rho;
      // Collect vectors along rays: start
      af::shared<double> vec;
      double p=0.;
      while(p<=1.) {
        double value = (int)(p * 100 + .5);
        p = (float)value / 100; // this make it precisely 0.00, 0.05, 0.10,...
        cctbx::cartesian<double> spc = cctbx::cartesian<>(
          center_cart[0] + p * (posc[0]-center_cart[0]),
          center_cart[1] + p * (posc[1]-center_cart[1]),
          center_cart[2] + p * (posc[2]-center_cart[2]) );
        cctbx::fractional<double> spf = unit_cell.fractionalize(spc);
        vec.push_back( tricubic_interpolation(map_data, spf) );
        p+=0.05;
      }
      vecs[i] = vec;
      // end
    }
    rho_mean = rho_mean/points_on_sphere_cart.size();
    rho_min_max_mean_ = af::tiny<double, 3>(rho_min, rho_max, rho_mean);
    //
    double cc_min=1.e+9, cc_max=-1.e+9, cc_mean=0.;
    int cntr=0;
    for(int i=0; i<vecs.size(); i++) {
      af::shared<double> vec_i = vecs[i];
      for(int j=0; j<vecs.size(); j++) {
        if(i<j) {
          af::shared<double> vec_j = vecs[j];
          double cc = scitbx::math::linear_correlation<>(
            vec_i.ref(), vec_j.ref()).coefficient();
          if(cc<cc_min) cc_min = cc;
          if(cc>cc_max) cc_max = cc;
          cc_mean += cc;
          cntr+=1;
        }
      }
    }
    cc_mean = cc_mean/cntr;
    ccs_min_max_mean_ = af::tiny<double, 3>(cc_min, cc_max, cc_mean);
    //
  }

  af::tiny<double, 3> rho_min_max_mean() {return rho_min_max_mean_;}
  af::tiny<double, 3> ccs_min_max_mean() {return ccs_min_max_mean_;}
};

class fit_point_3d_grid_search {
public:
  bool has_peak_;
  double map_best_, map_start_;
  cctbx::cartesian<> site_cart_moved_;
  fit_point_3d_grid_search(
    cctbx::cartesian<> const& site_cart,
    af::const_ref<double, af::c_grid<3> > const& map_data,
    double const& map_min, // TODO remove unused.
    cctbx::uctbx::unit_cell const& unit_cell,
    double const& amplitude,
    double const& increment)
  :
  has_peak_(true), site_cart_moved_(site_cart), map_best_(0), map_start_(0)
  {
    CCTBX_ASSERT(amplitude > 0.0 && increment > 0.0);
    double eps = 1.e-5;
    double x = site_cart[0];
    double y = site_cart[1];
    double z = site_cart[2];
    map_best_ = tricubic_interpolation(
      map_data, unit_cell.fractionalize(site_cart));
    map_start_ = map_best_;
    double x_shift = -amplitude;
    while(x_shift < amplitude) {
      x_shift += increment;
      double y_shift = -amplitude;
      double x_shifted = x+x_shift;
      while(y_shift < amplitude) {
        y_shift += increment;
        double y_shifted = y+y_shift;
        double z_shift = -amplitude;
        while(z_shift < amplitude) {
          z_shift += increment;
          cctbx::cartesian<> site_cart_ = cctbx::cartesian<>(
            x_shifted, y_shifted, z+z_shift);
          cctbx::fractional<> site_frac = unit_cell.fractionalize(site_cart_);
          double map_value = tricubic_interpolation(map_data, site_frac);
          if(map_value > map_best_) {
            map_best_ = map_value;
            site_cart_moved_ = site_cart_;
          }}}}
    double shift_x = std::abs(site_cart_moved_[0]-x);
    double shift_y = std::abs(site_cart_moved_[1]-y);
    double shift_z = std::abs(site_cart_moved_[2]-z);
    if(shift_x>amplitude || std::abs(shift_x-amplitude)<eps ||
       shift_y>amplitude || std::abs(shift_y-amplitude)<eps ||
       shift_z>amplitude || std::abs(shift_z-amplitude)<eps) {
      site_cart_moved_ = site_cart;
      has_peak_ = false;
    }
  }

  bool has_peak()                       {return has_peak_;}
  double map_best()                     {return map_best_;}
  double map_start()                    {return map_start_;}
  cctbx::cartesian<> site_cart_moved()  {return site_cart_moved_;}

};

template <typename DataType>
void
map_box_average(
  af::ref<DataType, af::c_grid<3> > map_data,
  DataType const& cutoff,
  int const& index_span)
{
  int nx = map_data.accessor()[0];
  int ny = map_data.accessor()[1];
  int nz = map_data.accessor()[2];
  for (int lx = 0; lx < nx; lx++) {
    for (int ly = 0; ly < ny; ly++) {
      for (int lz = 0; lz < nz; lz++) {
        if(map_data(lx,ly,lz)<cutoff) {
          DataType rho = 0.0;
          int counter = 0;
          for (int i = lx-index_span; i <= lx+index_span; i++) {
            int mx = i;
            if(i<0 || i>=nx) mx = scitbx::math::mod_positive(i, nx);
            for (int j = ly-index_span; j <= ly+index_span; j++) {
              int my = j;
              if(j<0 || j>=ny) my = scitbx::math::mod_positive(j, ny);
              for (int k = lz-index_span; k <= lz+index_span; k++) {
                int mz = k;
                if(k<0 || k>=nz) mz = scitbx::math::mod_positive(k, nz);
                rho += map_data(mx,my,mz);
                counter += 1;
          }}}
          map_data(lx,ly,lz) = rho / counter;
  }}}}
}

template <typename DataType>
af::shared<DataType> discrepancy_function(
  af::const_ref<DataType> const& map_1,
  af::const_ref<DataType> const& map_2,
  af::const_ref<DataType> const& cutoffs)
{
  CCTBX_ASSERT(af::max(map_1)<=1.);
  CCTBX_ASSERT(af::max(map_2)<=1.);
  CCTBX_ASSERT(af::min(map_1)>=0.);
  CCTBX_ASSERT(af::min(map_2)>=0.);
  CCTBX_ASSERT(af::min(cutoffs)>0. && af::max(cutoffs)<1.);
  CCTBX_ASSERT(map_1.size() == map_2.size());
  af::shared<DataType> result;
  int n_all = map_1.size();
  int n_diff = 0;
  for(int p = 0; p < cutoffs.size(); p++) {
    DataType q = cutoffs[p];
    for(int i = 0; i < map_1.size(); i++) {
      DataType rho1 = map_1[i];
      DataType rho2 = map_2[i];
      bool good = (rho1>=q&&rho2<q) || (rho1<q&&rho2>=q);
      if(good) n_diff += 1;
    }
    if(std::abs(1-q)>1.e-6 && std::abs(q)>1.e-6) {
      result.push_back(n_diff / (2*q*(1-q)*n_all));
    }
    n_diff = 0;
  }
  return result;
}

template <typename DataType>
af::shared<DataType> discrepancy_function(
  af::const_ref<DataType, af::c_grid<3> > const& map_1,
  af::const_ref<DataType, af::c_grid<3> > const& map_2,
  af::const_ref<DataType> const& cutoffs)
{
  CCTBX_ASSERT(af::max(map_1)<=1.);
  CCTBX_ASSERT(af::max(map_2)<=1.);
  CCTBX_ASSERT(af::min(map_1)>=0.);
  CCTBX_ASSERT(af::min(map_2)>=0.);
  CCTBX_ASSERT(af::min(cutoffs)>0. && af::max(cutoffs)<1.);
  af::c_grid<3> a1 = map_1.accessor();
  af::c_grid<3> a2 = map_2.accessor();
  for(int i = 0; i < 3; i++) CCTBX_ASSERT(a1[i]==a2[i]);
  af::shared<DataType> result;
  int n_all = a1[0]*a1[1]*a1[2];
  int n_diff = 0;
  for(int p = 0; p < cutoffs.size(); p++) {
    DataType q = cutoffs[p];
    for(int i = 0; i < a1[0]; i++) {
      for(int j = 0; j < a1[1]; j++) {
        for(int k = 0; k < a1[2]; k++) {
          DataType rho1 = map_1(i,j,k);
          DataType rho2 = map_2(i,j,k);
          bool good = (rho1>=q&&rho2<q) || (rho1<q&&rho2>=q);
          if(good) n_diff += 1;
        }
      }
    }
    if(std::abs(1-q)>1.e-6 && std::abs(q)>1.e-6) {
      result.push_back(n_diff / (2*q*(1-q)*n_all));
    }
    n_diff = 0;
  }
  return result;
}

template <typename DataType>
void
gamma_compression(
  af::ref<DataType, af::c_grid<3> > map_data,
  DataType const& gamma)
{
  CCTBX_ASSERT(gamma>0 && gamma<1);
  af::c_grid<3> a = map_data.accessor();
  for(int i = 0; i < a[0]; i++) {
    for(int j = 0; j < a[1]; j++) {
      for(int k = 0; k < a[2]; k++) {
         DataType rho = map_data(i,j,k);
         if(rho<0) map_data(i,j,k) = 0;
         else map_data(i,j,k) = std::pow(map_data(i,j,k), gamma);
  }}}
}

template <typename DataType>
void
sharpen(
  af::ref<DataType, af::c_grid<3> > map_data,
  int const& index_span,
  int const& n_averages,
  bool allow_negatives)
{
  af::c_grid<3> a = map_data.accessor();
  af::versa<double, af::c_grid<3> > result_map(a,
    af::init_functor_null<DataType>());
  af::ref<DataType, af::c_grid<3> > result_map_ref = result_map.ref();
  // copy
  for(int i = 0; i < a[0]; i++) {
    for(int j = 0; j < a[1]; j++) {
      for(int k = 0; k < a[2]; k++) {
         result_map_ref(i,j,k) = map_data(i,j,k);
  }}}
  // blur
  for(int i = 0; i < n_averages; i++) {
    map_box_average(result_map_ref, 9999., index_span);
  }
  //sharpen
  for(int i = 0; i < a[0]; i++) {
    for(int j = 0; j < a[1]; j++) {
      for(int k = 0; k < a[2]; k++) {
        if(allow_negatives) {
          map_data(i,j,k) = map_data(i,j,k)-result_map_ref(i,j,k);
        }
        else {
          map_data(i,j,k) = std::max(0., map_data(i,j,k)-result_map_ref(i,j,k));
        }
  }}}
}

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_UTILS_H
