#ifndef CCTBX_MAPTBX_UTILS_H
#define CCTBX_MAPTBX_UTILS_H

#include <cstddef>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <cctbx/uctbx.h>
#include <scitbx/math/utils.h>
#include <scitbx/math/modulo.h>
#include <scitbx/array_family/sort.h>

#include <cctbx/maptbx/eight_point_interpolation.h> //indirect import?

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

template <typename DataType>
void set_box(
  af::const_ref<DataType, af::c_grid<3> > const& map_data_from,
  af::ref<DataType, af::c_grid<3> > map_data_to,
  af::tiny<int, 3> const& start,
  af::tiny<int, 3> const& end)
{
  int ii=0;
  for (int i = start[0]; i < end[0]; i++) {
    int jj=0;
    for (int j = start[1]; j < end[1]; j++) {
      int kk=0;
      for (int k = start[2]; k < end[2]; k++) {
        map_data_to(i,j,k) = map_data_from(ii,jj,kk);
        kk+=1;
      }
      jj+=1;
    }
    ii+=1;
  }
}

template <typename DataType>
void cut_by(
       af::ref<DataType, af::c_grid<3> > kick,
       af::ref<DataType, af::c_grid<3> > fem,
       DataType cut_by_threshold)
{
  af::tiny<int, 3> a1 = kick.accessor();
  af::tiny<int, 3> a2 = fem.accessor();
  for(int i = 0; i < 3; i++) CCTBX_ASSERT(a1[i]==a2[i]);
  for(int i = 0; i < a1[0]; i++) {
    for(int j = 0; j < a1[1]; j++) {
      for(int k = 0; k < a1[2]; k++) {
         if(kick(i,j,k)<cut_by_threshold) {
           kick(i,j,k)=0;
  }}}}
  for(int i = 0; i < a1[0]; i++) {
    for(int j = 0; j < a1[1]; j++) {
      for(int k = 0; k < a1[2]; k++) {
         if(fem(i,j,k)<1) {
           fem(i,j,k)=0;
  }}}}

  for(int i = 0; i < a1[0]; i++) {
    for(int j = 0; j < a1[1]; j++) {
      for(int k = 0; k < a1[2]; k++) {
         double rk = kick(i,j,k);
         double rf = fem(i,j,k);
         if(rk==0 || rf==0) {
           kick(i,j,k)=0;
           fem(i,j,k)=0;
  }}}}
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
       DataType threshold)
{
  af::tiny<int, 3> a1 = map_data_1.accessor();
  af::tiny<int, 3> a2 = map_data_2.accessor();
  for(int i = 0; i < 3; i++) CCTBX_ASSERT(a1[i]==a2[i]);
  for(int i = 0; i < a1[0]; i++) {
    for(int j = 0; j < a1[1]; j++) {
      for(int k = 0; k < a1[2]; k++) {
         double rho1 = map_data_1(i,j,k);
         double rho2 = map_data_2(i,j,k);
         bool c1 = rho1>threshold && rho2<threshold;
         bool c2 = rho2>threshold && rho1<threshold;
         if(c1 || c2) {
           map_data_1(i,j,k)=0;
           map_data_2(i,j,k)=0;
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
  DataType const& cutoff,
  int const& index_span)
{
  int nx = map_data.accessor()[0];
  int ny = map_data.accessor()[1];
  int nz = map_data.accessor()[2];
  for (int lx = 0; lx < nx; lx++) {
    for (int ly = 0; ly < ny; ly++) {
      for (int lz = 0; lz < nz; lz++) {
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
            for (int j = ly-index_span; j <= ly+index_span; j++) {
              for (int k = lz-index_span; k <= lz+index_span; k++) {
                int mx = scitbx::math::mod_positive(i, nx);
                int my = scitbx::math::mod_positive(j, ny);
                int mz = scitbx::math::mod_positive(k, nz);
                rho += map_data(mx,my,mz);
                counter += 1;
          }}}
          map_data(lx,ly,lz) = rho / counter;
  }}}}
}

template <typename DataType>
void
sharpen(
  af::ref<DataType, af::c_grid<3> > map_data,
  int const& index_span,
  int const& n_averages)
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
        map_data(i,j,k) = std::max(0., map_data(i,j,k)-result_map_ref(i,j,k));
  }}}
}

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_UTILS_H
