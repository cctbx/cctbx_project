#ifndef CCTBX_MAPTBX_AVERAGE_DENSITIES_H
#define CCTBX_MAPTBX_AVERAGE_DENSITIES_H

#include <cctbx/uctbx.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/math/fast_approx_math.h>
#include <scitbx/math/utils.h>
#include <cctbx/maptbx/interpolation.h>
#include <scitbx/math/linear_correlation.h>
#include <cctbx/maptbx/histogram.h>

namespace cctbx { namespace maptbx {

af::versa<double, af::c_grid<3> > denmod_simple(
  af::const_ref<double, af::c_grid<3> > const& map_data,
  af::tiny<int, 3> const& n_real,
  double cutoffp,
  double cutoffm)
{
  int nx = n_real[0];
  int ny = n_real[1];
  int nz = n_real[2];
  af::versa<double, af::c_grid<3> > result_map(af::c_grid<3>(nx,ny,nz),
    af::init_functor_null<double>());
  af::ref<double, af::c_grid<3> > result_map_ref = result_map.ref();
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        double m = map_data(i,j,k);
        if(m > cutoffp) {
          result_map_ref(i,j,k) = m-cutoffp;
        }
        else if(m < cutoffm) {
          result_map_ref(i,j,k) = -m+cutoffm;
        }
        else {
          result_map_ref(i,j,k) = 0;
        }
        CCTBX_ASSERT(result_map_ref(i,j,k) >= 0);
  }}}
  return result_map;
}

af::versa<double, af::c_grid<3> > combine_and_maximize_maps(
                   af::const_ref<double, af::c_grid<3> > const& map_data_1,
                   af::const_ref<double, af::c_grid<3> > const& map_data_2,
                   af::tiny<int, 3> const& n_real)
{
  int nx = n_real[0];
  int ny = n_real[1];
  int nz = n_real[2];
  af::versa<double, af::c_grid<3> > result_map(af::c_grid<3>(nx,ny,nz),
    af::init_functor_null<double>());
  af::ref<double, af::c_grid<3> > result_map_ref = result_map.ref();
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        double m1 = map_data_1(i,j,k);
        double m2 = map_data_2(i,j,k);
        if(std::abs(m1) >= std::abs(m2)) result_map_ref(i,j,k) = m1;
        else result_map_ref(i,j,k) = m2;
  }}}
  return result_map;
}

af::versa<double, af::c_grid<3> > superpose_maps(
                   cctbx::uctbx::unit_cell const& unit_cell_1,
                   cctbx::uctbx::unit_cell const& unit_cell_2,
                   af::const_ref<double, af::c_grid<3> > const& map_data_1,
                   af::tiny<int, 3> const& n_real_2,
                   scitbx::mat3<double> const& rotation_matrix,
                   scitbx::vec3<double> const& translation_vector,
                   bool wrapping)
{
  int nx = n_real_2[0];
  int ny = n_real_2[1];
  int nz = n_real_2[2];
  af::versa<double, af::c_grid<3> > result_map(af::c_grid<3>(nx,ny,nz),
    af::init_functor_null<double>());
  af::ref<double, af::c_grid<3> > result_map_ref = result_map.ref();
  for (int i = 0; i < nx; i++) {
    double gpfx = i/static_cast<double>(nx);
    for (int j = 0; j < ny; j++) {
      double gpfy = j/static_cast<double>(ny);
      double lb = -1.e-6;
      double ub = 1+1.e-6;
      for (int k = 0; k < nz; k++) {
        cctbx::fractional<> grid_frac_in_1 = cctbx::fractional<>(gpfx, gpfy,
          k/static_cast<double>(nz));
        cctbx::cartesian<> grid_cart_in_1 = unit_cell_2.orthogonalize(
          grid_frac_in_1);
        cctbx::cartesian<> point_cart_in_2 = rotation_matrix * grid_cart_in_1 +
          translation_vector;
        cctbx::fractional<> point_frac_in_2 =
          unit_cell_1.fractionalize(point_cart_in_2);
        if (wrapping){
          result_map_ref(i,j,k) = tricubic_interpolation(map_data_1,
            point_frac_in_2);
        } else { // Separate cases for inside (0,1), very very close, further
          if (
             point_frac_in_2[0] >= 0 && point_frac_in_2[0] <= 1 &&
             point_frac_in_2[1] >= 0 && point_frac_in_2[1] <= 1 &&
             point_frac_in_2[2] >= 0 && point_frac_in_2[2] <= 1 ) {
            // cleanly inside
            result_map_ref(i,j,k) = tricubic_interpolation(map_data_1,
              point_frac_in_2);
          } else if (
             point_frac_in_2[0]<lb || point_frac_in_2[0] > ub ||
             point_frac_in_2[1]<lb || point_frac_in_2[1] > ub ||
             point_frac_in_2[2]<lb || point_frac_in_2[2] > ub ) {
            // cleanly outside (further than 1.e-6 outside boundaries of (0,1)
            result_map_ref(i,j,k) = 0.;
          } else {
            // just outside. Move to exact boundary
            for (int i = 0; i < 3; i++) {
              if (point_frac_in_2[i]>lb && point_frac_in_2[i]< 0 ){
                point_frac_in_2[i] = 0;}
              if (point_frac_in_2[i]>1 && point_frac_in_2[i]< ub ){
                point_frac_in_2[i] = 1;}
            }
            result_map_ref(i,j,k) = tricubic_interpolation(map_data_1,
              point_frac_in_2);
          }
        }
  }}}
  return result_map;
}

af::versa<double, af::c_grid<3> > rotate_translate_map(
                   cctbx::uctbx::unit_cell const& unit_cell,
                   af::const_ref<double, af::c_grid<3> > const& map_data,
                   scitbx::mat3<double> const& rotation_matrix,
                   scitbx::vec3<double> const& translation_vector,
                   af::tiny<int, 3> const& start,
                   af::tiny<int, 3> const& end)
{
    int nx = static_cast<int>(map_data.accessor()[0]);
    int ny = static_cast<int>(map_data.accessor()[1]);
    int nz = static_cast<int>(map_data.accessor()[2]);
    af::versa<double, af::c_grid<3> > new_data(af::c_grid<3>(nx,ny,nz),
      af::init_functor_null<double>());
    af::ref<double, af::c_grid<3> > new_data_ref = new_data.ref();
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
          if(i>=start[0] && j>=start[1] && k>=start[2] &&
             i<=end[0] && j<=end[1] && k<=end[2]) {
            cctbx::fractional<> grid_node_frac = cctbx::fractional<>(
              i/static_cast<double>(nx),
              j/static_cast<double>(ny),
              k/static_cast<double>(nz));
            cctbx::cartesian<> grid_node_cart = unit_cell.orthogonalize(
              grid_node_frac);
            cctbx::fractional<> grid_node_frac_shifted = unit_cell.fractionalize(
              rotation_matrix * grid_node_cart + translation_vector);
            for (int p = 0; p < 5; p++) {
              for (int q = 0; q < 3; q++) {
                if(grid_node_frac_shifted[q] <  0) grid_node_frac_shifted[q] += 1;
                if(grid_node_frac_shifted[q] >= 1) grid_node_frac_shifted[q] -= 1;
              }
            }
            new_data_ref(i,j,k) = tricubic_interpolation(map_data,
              grid_node_frac_shifted);
          }
    }}}
    return new_data;
}

af::versa<double, af::c_grid<3> > rotate_translate_map(
                   cctbx::uctbx::unit_cell const& unit_cell,
                   af::const_ref<double, af::c_grid<3> > const& map_data,
                   scitbx::mat3<double> const& rotation_matrix,
                   scitbx::vec3<double> const& translation_vector)
{
    int nx = static_cast<int>(map_data.accessor()[0]);
    int ny = static_cast<int>(map_data.accessor()[1]);
    int nz = static_cast<int>(map_data.accessor()[2]);
    af::versa<double, af::c_grid<3> > new_data(af::c_grid<3>(nx,ny,nz),
      af::init_functor_null<double>());
    af::ref<double, af::c_grid<3> > new_data_ref = new_data.ref();
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
          cctbx::fractional<> grid_node_frac = cctbx::fractional<>(
            i/static_cast<double>(nx),
            j/static_cast<double>(ny),
            k/static_cast<double>(nz));
          cctbx::cartesian<> grid_node_cart = unit_cell.orthogonalize(
            grid_node_frac);
          cctbx::fractional<> grid_node_frac_shifted = unit_cell.fractionalize(
            rotation_matrix * grid_node_cart + translation_vector);
          for (int p = 0; p < 5; p++) {
            for (int q = 0; q < 3; q++) {
              if(grid_node_frac_shifted[q] <  0) grid_node_frac_shifted[q] += 1;
              if(grid_node_frac_shifted[q] >= 1) grid_node_frac_shifted[q] -= 1;
            }
          }
          new_data_ref(i,j,k) = tricubic_interpolation(map_data,
            grid_node_frac_shifted);
    }}}
    return new_data;
}

void
  cc_weighted_maps(
    af::ref<double, af::c_grid<3> > data_1,
    af::ref<double, af::c_grid<3> > data_2)
{
  int NX = data_1.accessor()[0];
  int NY = data_1.accessor()[1];
  int NZ = data_1.accessor()[2];
  // perhaps need pre-filter maps first to avoid computing CC(0,0)
  for (int lx = 0; lx < NX; lx++) {
    for (int ly = 0; ly < NY; ly++) {
      for (int lz = 0; lz < NZ; lz++) {
        af::shared<double> m1;
        af::shared<double> m2;
        int inc = 1;
        for (int i = lx-inc; i <= lx+inc; i++) {
          for (int j = ly-inc; j <= ly+inc; j++) {
            for (int k = lz-inc; k <= lz+inc; k++) {
              int mx = scitbx::math::mod_positive(i, NX);
              int my = scitbx::math::mod_positive(j, NY);
              int mz = scitbx::math::mod_positive(k, NZ);
              m1.push_back(data_1(mx,my,mz));
              m2.push_back(data_2(mx,my,mz));
        }}}
        double cc =
          scitbx::math::linear_correlation<>(m1.ref(), m2.ref()).coefficient();
        data_1(lx,ly,lz) = data_1(lx,ly,lz)*cc;
        data_2(lx,ly,lz) = data_2(lx,ly,lz)*cc;
  }}}
}

class volume_scale {
public:
  af::versa<double, af::c_grid<3> > map_new;
  af::shared<double> v_values_;

  volume_scale(
    af::const_ref<double, af::c_grid<3> > const& map,
    int const& n_bins)
  {
    int nx = map.accessor()[0];
    int ny = map.accessor()[1];
    int nz = map.accessor()[2];
    map_new.resize(af::c_grid<3>(nx,ny,nz), 0);
    double rho_min = af::min(map);
    histogram hist = histogram(map, n_bins);
    double bin_width = hist.bin_width();
    v_values_ = hist.c_values();
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
          double rho = map(i,j,k);
          int index = scitbx::math::nearest_integer((rho-rho_min)/bin_width);
          if(index<0) index=0;
          if(index>=n_bins) index=n_bins-1;
          double rho_new = 0;
          if(index+1<n_bins) {
            double rho_n = rho_min + index*bin_width;
            rho_new = v_values_[index] +
              (v_values_[index+1]-v_values_[index]) * (rho-rho_n)/bin_width;
            if(rho_new<0) rho_new = v_values_[index];
          }
          else {
            rho_new = v_values_[index];
          }
          CCTBX_ASSERT(rho_new>=0);
          map_new(i,j,k) = rho_new;
        }
      }
    }
  }

  af::versa<double, af::c_grid<3> > map_data() {return map_new;}
  af::shared<double> v_values() {return v_values_;}

};

class volume_scale_1d {
public:
  af::shared<double> map_new;
  af::shared<double> v_values_;

  volume_scale_1d(
    af::const_ref<double> const& map,
    int const& n_bins)
  {
    map_new.resize(map.size());
    map_new.fill(0);
    double rho_min = af::min(map);
    histogram hist = histogram(map, n_bins);
    double bin_width = hist.bin_width();
    v_values_ = hist.c_values();
    for(int i = 0; i < map.size(); i++) {
      double rho = map[i];
      int index = scitbx::math::nearest_integer((rho-rho_min)/bin_width);
      if(index<0) index=0;
      if(index>=n_bins) index=n_bins-1;
      double rho_new = 0;
      if(index+1<n_bins) {
        double rho_n = rho_min + index*bin_width;
        rho_new = v_values_[index] +
          (v_values_[index+1]-v_values_[index]) * (rho-rho_n)/bin_width;
        if(rho_new<0) rho_new = v_values_[index];
      }
      else {
        rho_new = v_values_[index];
      }
      CCTBX_ASSERT(rho_new>=0);
      map_new[i] = rho_new;
    }
  }

  af::shared<double> map_data() {return map_new;}
  af::shared<double> v_values() {return v_values_;}

};

class non_linear_map_modification_to_match_average_cumulative_histogram {
public:
  af::versa<double, af::c_grid<3> > map_1_new;
  af::versa<double, af::c_grid<3> > map_2_new;
  af::shared<double> sigs;
  af::shared<double> h1;
  af::shared<double> h2;
  af::shared<double> h12;

  non_linear_map_modification_to_match_average_cumulative_histogram(
    af::const_ref<double, af::c_grid<3> > const& map_1,
    af::const_ref<double, af::c_grid<3> > const& map_2)
  {
    int nx1 = map_1.accessor()[0];
    int ny1 = map_1.accessor()[1];
    int nz1 = map_1.accessor()[2];
    int nx2 = map_2.accessor()[0];
    int ny2 = map_2.accessor()[1];
    int nz2 = map_2.accessor()[2];
    CCTBX_ASSERT(nx1==nx2 && ny1==ny2 && nz1==nz2);
    map_1_new.resize(af::c_grid<3>(nx1,ny1,nz1), 0);
    map_2_new.resize(af::c_grid<3>(nx1,ny1,nz1), 0);
    double max1 = af::max(map_1);
    double min1 = af::min(map_1);
    double max2 = af::max(map_2);
    double min2 = af::min(map_2);
    double mean1 = 0;
    double mean2 = 0;
    for (int i = 0; i < nx1; i++) {
      for (int j = 0; j < ny1; j++) {
        for (int k = 0; k < nz1; k++) {
          // do not shift to be between 0 and 1: operate on map as
          //double m2 = (map_2(i,j,k)-min2)*(2)/(max2-min2)-1;
          //double m1 = (map_1(i,j,k)-min1)*(2)/(max1-min1)-1;
          //map_2_new(i,j,k)=m2;
          //map_1_new(i,j,k)=m1;
          //mean1 = mean1 + m1;
          //mean2 = mean2 + m2;
          map_2_new(i,j,k)=map_2(i,j,k);
          map_1_new(i,j,k)=map_1(i,j,k);
        }
      }
    }
    // do not subtract mean
    //CCTBX_ASSERT(std::abs(af::max(map_1_new.ref())-af::max(map_2_new.ref()))<1e-6);
    //CCTBX_ASSERT(std::abs(af::min(map_1_new.ref())-af::min(map_2_new.ref()))<1e-6);
    //mean1 = mean1/map_1.size();
    //mean2 = mean2/map_1.size();
    //for (int i = 0; i < nx1; i++) {
    //  for (int j = 0; j < ny1; j++) {
    //    for (int k = 0; k < nz1; k++) {
    //      map_2_new(i,j,k)=map_2_new(i,j,k) - mean2;
    //      map_1_new(i,j,k)=map_1_new(i,j,k) - mean1;
    //    }
    //  }
    //}
    int n_bins = 3000;
    max1 = af::max(map_1_new.ref());
    min1 = af::min(map_1_new.ref());
    max2 = af::max(map_2_new.ref());
    min2 = af::min(map_2_new.ref());
    double start = std::min(min1,min2);
    double end   = std::max(max1,max2);
    histogram hist1 = histogram(map_1_new.ref(), n_bins, start, end);
    histogram hist2 = histogram(map_2_new.ref(), n_bins, start, end);
    h1 = hist1.c_values();
    h2 = hist2.c_values();
    sigs = hist1.arguments(); // are sigs from hist2 the same?
    double bin_width = hist1.bin_width();
    // can use mean of two or one of the two:
    //for (int i = 0; i < h1.size(); i++) {
    //  h12.push_back((h1[i]+h2[i])/2);
    //}
    for (int i = 0; i < h1.size(); i++) {
      h12.push_back(h2[i]);
    }

    map_mod(map_1_new.ref(), map_2_new.ref(), h1.ref(), h2.ref(), sigs.ref(),
      h12.ref(), bin_width, start);
  }

  af::versa<double, af::c_grid<3> > map_1() {return map_1_new;}
  af::versa<double, af::c_grid<3> > map_2() {return map_2_new;}
  af::shared<double> histogram_1() {return h1;}
  af::shared<double> histogram_2() {return h2;}
  af::shared<double> histogram_12() {return h12;}
  af::shared<double> histogram_values() {return sigs;}

protected:
  void map_mod(
    af::ref<double, af::c_grid<3> > m1,
    af::ref<double, af::c_grid<3> > m2,
    af::const_ref<double> h1,
    af::const_ref<double> h2,
    af::const_ref<double> sigs,
    af::const_ref<double> h_ave,
    double bw,
    double start) {
    int nx = m1.accessor()[0];
    int ny = m1.accessor()[1];
    int nz = m1.accessor()[2];
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
          double h1_best = h1[scitbx::math::nearest_integer((m1(i,j,k)-start)/bw)];
          double h2_best = h2[scitbx::math::nearest_integer((m2(i,j,k)-start)/bw)];
          //
          double sigav1_best = sigs[0];
          double hav1_best = 0;
          double r1 = 1e+9;
          double sigav2_best = sigs[0];
          double hav2_best = 0;
          double r2 = 1e+9;
          for (int n=0; n<h1.size(); n++) {
            double r1_ = std::abs(h1_best-h_ave[n]);
            double r2_ = std::abs(h2_best-h_ave[n]);
            if(r1_<r1) {
              r1 = r1_;
              sigav1_best = sigs[n];
              hav1_best = h_ave[n];
            }
            if(r2_<r2) {
              r2 = r2_;
              sigav2_best = sigs[n];
              hav2_best = h_ave[n];
            }
          }
          m1(i,j,k)=sigav1_best;
          m2(i,j,k)=sigav2_best;
        }
      }
    }
  }
};

class grid_points_in_sphere_around_atom_and_distances {
public:
  grid_points_in_sphere_around_atom_and_distances(
                           cctbx::uctbx::unit_cell const& uc,
                           af::const_ref<double, af::c_grid<3> > const& data,
                           double const& rad,
                           double const& shell,
                           scitbx::vec3<double> const& site_frac)
  {
    int nx = data.accessor()[0];
    int ny = data.accessor()[1];
    int nz = data.accessor()[2];
    double distsq;
    double grid_step = uc.parameters()[3] / nx;
    int mx,my,mz,kx,ky,kz;
    double abc = uc.parameters()[0]*uc.parameters()[1]*uc.parameters()[2];
    double uc_volume_over_abc = uc.volume() / abc;
    double sin_alpha = std::sin(scitbx::deg_as_rad(uc.parameters()[3]));
    double sin_beta  = std::sin(scitbx::deg_as_rad(uc.parameters()[4]));
    double sin_gamma = std::sin(scitbx::deg_as_rad(uc.parameters()[5]));
    double ascale = uc_volume_over_abc / sin_alpha;
    double bscale = uc_volume_over_abc / sin_beta;
    double cscale = uc_volume_over_abc / sin_gamma;
    double as = uc.parameters()[0]/ascale;
    double bs = uc.parameters()[1]/bscale;
    double cs = uc.parameters()[2]/cscale;
    double xshell = shell/uc.parameters()[0]/ascale;
    double yshell = shell/uc.parameters()[1]/bscale;
    double zshell = shell/uc.parameters()[2]/cscale;
    double mr1= uc.metrical_matrix()[0]; //a*a;
    double mr2= uc.metrical_matrix()[3]; //a*b*cos(gamma)
    double mr3= uc.metrical_matrix()[4]; //a*c*cos(beta)
    double mr5= uc.metrical_matrix()[1]; //b*b
    double mr6= uc.metrical_matrix()[5]; //c*b*cos(alpha)
    double mr9= uc.metrical_matrix()[2]; //c*c
    double tmr2 = mr2*2.0; //2*a*b*cos(gamma);
    double tmr3 = mr3*2.0; //2*a*c*cos(beta);
    double tmr6 = mr6*2.0; //2*b*c*cos(alpha);
    double xL = -xshell;
    double xR = 1.0+xshell;
    double yL = -yshell;
    double yR = 1.0+yshell;
    double zL = -zshell;
    double zR = 1.0+zshell;
    double xf = site_frac[0];
    double yf = site_frac[1];
    double zf = site_frac[2];
    double radsq = rad * rad;
    bool close_to_au=(xf>=xL||xf<=xR)&&(yf>=yL||yf<=yR)&&(zf>=zL||zf<=zR);
    if(close_to_au) {
       double coas = rad/as;
       double cobs = rad/bs;
       double cocs = rad/cs;
       int x1box = scitbx::math::nearest_integer( nx*(xf-coas) ) - 1;
       int x2box = scitbx::math::nearest_integer( nx*(xf+coas) ) + 1;
       int y1box = scitbx::math::nearest_integer( ny*(yf-cobs) ) - 1;
       int y2box = scitbx::math::nearest_integer( ny*(yf+cobs) ) + 1;
       int z1box = scitbx::math::nearest_integer( nz*(zf-cocs) ) - 1;
       int z2box = scitbx::math::nearest_integer( nz*(zf+cocs) ) + 1;
       for(kx = x1box; kx <= x2box; kx++) {
           double xn=xf-double(kx)/nx;
           for(ky = y1box; ky <= y2box; ky++) {
               double yn=yf-double(ky)/ny;
               for(kz = z1box; kz <= z2box; kz++) {
                   double zn=zf-double(kz)/nz;
                   distsq=mr1*xn*xn+mr5*yn*yn+mr9*zn*zn+
                                              tmr2*xn*yn+tmr3*xn*zn+tmr6*yn*zn;
                   if(distsq <= radsq) {
                      mx = scitbx::math::mod_positive(kx, nx);
                      my = scitbx::math::mod_positive(ky, ny);
                      mz = scitbx::math::mod_positive(kz, nz);
                      //double rho = 2.0 * std::exp (- 3.0 * distsq );
                      data_at_grid_points_.push_back(data(mx,my,mz));
                      //data_at_grid_points_.push_back(rho);
                      distances_.push_back(std::sqrt(distsq));
                   }
               }
           }
       }
    }
    //
    double tolerance = grid_step/25. ;
    for(std::size_t i = 0; i < data_at_grid_points_.size(); i++) {
        double dist = distances_[i];
        double data_ave = data_at_grid_points_[i];
        int counter = 1;
        for(std::size_t j = 0; j < data_at_grid_points_.size(); j++) {
            if(distances_[j]<dist+tolerance && distances_[j]>dist-tolerance && i!=j && std::abs(distances_[i]-distances_[j]) > 1.e-6) {
               counter++;
               data_ave = data_ave + data_at_grid_points_[j];
//               if (distances_[j] < 0.00001)
//std::cout<<data_at_grid_points_[i]<<" "<<data_at_grid_points_[j]<<" "<<counter<<" "<<data_ave/counter<<" "<<distances_[j]<<" "<<i<<" "<<j<<dist<<std::endl;
            }
        }
        data_at_grid_points_averaged_.push_back(data_ave / counter);
    }
    //
  }
  af::shared<double> data_at_grid_points() { return data_at_grid_points_; }
  af::shared<double> data_at_grid_points_averaged() { return data_at_grid_points_averaged_; }
  af::shared<double> distances() { return distances_; }
protected:
  af::shared<double> data_at_grid_points_, data_at_grid_points_averaged_;
  af::shared<double> distances_;
};

class find_gaussian_parameters {
public:
  find_gaussian_parameters() {}
  find_gaussian_parameters(af::const_ref<double> const& data_at_grid_points,
                           af::const_ref<double> const& distances,
                           double const& cutoff_radius,
                           double const& allowed_region_radius,
                           double weight_power)
  {
    CCTBX_ASSERT(cutoff_radius <= allowed_region_radius);
    double max_distances = af::max(distances);
    CCTBX_ASSERT(cutoff_radius <= max_distances &&
                 allowed_region_radius <= max_distances);
    double p=0.0, q=0.0, r=0.0, s=0.0, t=0.0;
    int n = 0;
    int number_of_points = data_at_grid_points.size();
    for(int i = 0; i < number_of_points; i++) {
        if(data_at_grid_points[i] > 0.0 && distances[i] <= cutoff_radius) {
           n = n + 1;
           double data_i = data_at_grid_points[i];
           double distance_i = distances[i];
           double distance_i_sq = distance_i*distance_i;
           double pow_distance_i = std::pow(distance_i, weight_power);
           double log_of_data_i = std::log(data_i);
           if(pow_distance_i < 1.e-9) pow_distance_i = 1.0;
           p=p+log_of_data_i/pow_distance_i;
           q=q+distance_i_sq/pow_distance_i;
           r=r+(distance_i_sq*distance_i_sq)/pow_distance_i;
           s=s+distance_i_sq*log_of_data_i/pow_distance_i;
           t=t+1.0/pow_distance_i;
        }
    }
    CCTBX_ASSERT(r != 0.0);
    double den = t-q*q/r;
    CCTBX_ASSERT(den != 0.0);
    double alpha_opt = (p-s*q/r) / den;
    a_real_space_ = std::exp( alpha_opt );
    b_real_space_ = 1./r*(alpha_opt*q - s);
    double tmp = b_real_space_/scitbx::constants::pi;
    double pi_sq = scitbx::constants::pi*scitbx::constants::pi;
    CCTBX_ASSERT(tmp != 0.0);
    a_reciprocal_space_ = a_real_space_/std::sqrt(tmp*tmp*tmp);
    CCTBX_ASSERT(b_real_space_ != 0.0);
    b_reciprocal_space_ = pi_sq/b_real_space_*4; // not deja divise par 4 !!!
    double num = 0.0;
    double denum = 0.0;
    for(int i = 0; i < number_of_points; i++) {
        if(data_at_grid_points[i]>0.0 && distances[i]<=allowed_region_radius) {
           double data_i = data_at_grid_points[i];
           double distance_i = distances[i];
           double distance_i_sq = distance_i*distance_i;
           double data_i_approx =
                        a_real_space_ * std::exp(-b_real_space_*distance_i_sq);
            num=num+std::abs(data_i-data_i_approx);
            denum=denum+data_i;
        }
    }
    CCTBX_ASSERT(denum != 0.0);
    gof_ = num/denum*100.;
  }
  double a_real_space() { return a_real_space_; }
  double b_real_space() { return b_real_space_; }
  double a_reciprocal_space() { return a_reciprocal_space_; }
  double b_reciprocal_space() { return b_reciprocal_space_; }
  double gof() { return gof_; }
protected:
  double a_real_space_,b_real_space_,a_reciprocal_space_,b_reciprocal_space_;
  double gof_;
};

class one_gaussian_peak_approximation {
public:
  one_gaussian_peak_approximation(
                              af::const_ref<double> const& data_at_grid_points,
                              af::const_ref<double> const& distances,
                              bool const& use_weights,
                              bool const& optimize_cutoff_radius)
  {
    first_zero_radius_ =find_first_zero_radius(data_at_grid_points, distances);
    double radius_increment = 0.01, weight_power = 0.0, weight_increment =0.05;
    weight_power_ = -1;
    cutoff_radius_ = -1;
    gof_ = 999.;
    if(use_weights && optimize_cutoff_radius) {
       for(double c_radius = 0.1; c_radius <= first_zero_radius_;
                                        c_radius = c_radius+radius_increment) {
       for(double w_power=0.0;w_power<=10.0;w_power=w_power+weight_increment) {
           find_gaussian_parameters fgp_obj(data_at_grid_points,
                                            distances,
                                            c_radius,
                                            first_zero_radius_,
                                            w_power);
           if(fgp_obj.gof() < gof_) {
              gof_           = fgp_obj.gof();
              weight_power_  = w_power;
              cutoff_radius_ = c_radius;
              fgp_obj_       = fgp_obj;
           }
       }}

    }
    else if(use_weights && !optimize_cutoff_radius) {
       for(double w_power=0.0;w_power<=20.0;w_power=w_power+weight_increment) {
           find_gaussian_parameters fgp_obj(data_at_grid_points,
                                            distances,
                                            first_zero_radius_,
                                            first_zero_radius_,
                                            w_power);
           if(fgp_obj.gof() < gof_) {
              gof_           = fgp_obj.gof();
              weight_power_  = w_power;
              cutoff_radius_ = first_zero_radius_;
              fgp_obj_       = fgp_obj;
           }
       }
    }
    else if(optimize_cutoff_radius && !use_weights) {
       for(double c_radius = 0.1; c_radius <= first_zero_radius_;
                                        c_radius = c_radius+radius_increment) {
           find_gaussian_parameters fgp_obj(data_at_grid_points,
                                            distances,
                                            c_radius,
                                            first_zero_radius_,
                                            0.0);
           if(fgp_obj.gof() < gof_) {
              gof_           = fgp_obj.gof();
              cutoff_radius_ = c_radius;
              fgp_obj_       = fgp_obj;
           }
       }
    }
    else {
       find_gaussian_parameters fgp_obj(data_at_grid_points,
                                        distances,
                                        first_zero_radius_,
                                        first_zero_radius_,
                                        0.0);
       gof_           = fgp_obj.gof();
       cutoff_radius_ = first_zero_radius_;
       fgp_obj_       = fgp_obj;
    }
  }
  double find_first_zero_radius(
                              af::const_ref<double> const& data_at_grid_points,
                              af::const_ref<double> const& distances)
  {
    double shift = 1.e-3;
    af::shared<double> distances_for_negative_data_at_grid_points;
    for(int i = 0; i < data_at_grid_points.size(); i++) {
      if(data_at_grid_points[i] < 0.0) {
        distances_for_negative_data_at_grid_points.push_back(distances[i]);
      }
    }
    double result = af::max(distances);
    if(distances_for_negative_data_at_grid_points.size() > 0) {
       result =
         af::min(distances_for_negative_data_at_grid_points.const_ref())-shift;
    }
    return result;
  }
  double a_real_space()       { return fgp_obj_.a_real_space(); }
  double b_real_space()       { return fgp_obj_.b_real_space(); }
  double a_reciprocal_space() { return fgp_obj_.a_reciprocal_space(); }
  double b_reciprocal_space() { return fgp_obj_.b_reciprocal_space(); }
  double gof()
  {
   CCTBX_ASSERT(gof_ == fgp_obj_.gof());
   return gof_;
  }
  double cutoff_radius()      { return cutoff_radius_; }
  double weight_power()       { return weight_power_; }
  double first_zero_radius()  { return first_zero_radius_; }
protected:
  double a_real_space_,b_real_space_,a_reciprocal_space_,b_reciprocal_space_;
  double gof_, cutoff_radius_, weight_power_, first_zero_radius_;
  find_gaussian_parameters fgp_obj_;
};

  template <typename DataType>
  af::shared<DataType>
  average_densities(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<DataType, af::c_grid<3> > const& data,
    af::const_ref<scitbx::vec3<double> > const& sites_frac,
    float radius)
  {
    // Severe code duplication: mmtbx/masks/around_atoms.h
    af::shared<DataType> result((af::reserve(sites_frac.size())));
    typedef float f_t;
    int nx = static_cast<int>(data.accessor()[0]);
    int ny = static_cast<int>(data.accessor()[1]);
    int nz = static_cast<int>(data.accessor()[2]);
    f_t mr1= static_cast<f_t>(unit_cell.metrical_matrix()[0]); // a*a
    f_t mr5= static_cast<f_t>(unit_cell.metrical_matrix()[1]); // b*b
    f_t mr9= static_cast<f_t>(unit_cell.metrical_matrix()[2]); // c*c
    // a*b*cos(gamma)
    f_t mr2= static_cast<f_t>(unit_cell.metrical_matrix()[3]);
    // a*c*cos(beta)
    f_t mr3= static_cast<f_t>(unit_cell.metrical_matrix()[4]);
    // c*b*cos(alpha)
    f_t mr6= static_cast<f_t>(unit_cell.metrical_matrix()[5]);
    f_t tmr2 = mr2*2; //2*a*b*cos(gamma);
    f_t tmr3 = mr3*2; //2*a*c*cos(beta);
    f_t tmr6 = mr6*2; //2*b*c*cos(alpha);
    f_t sx = 1/static_cast<f_t>(nx); f_t tsx= sx*2; f_t sxsq=mr1*sx*sx;
    f_t sy = 1/static_cast<f_t>(ny); f_t tsy= sy*2; f_t sysq=mr5*sy*sy;
    f_t sz = 1/static_cast<f_t>(nz); f_t tsz= sz*2; f_t szsq=mr9*sz*sz;
    f_t w1=mr1*sx*tsx; f_t w4=mr5*sy*tsy;
    f_t w2=mr2*sx*tsy; f_t w5=mr6*sy*tsz;
    f_t w3=mr3*sx*tsz; f_t w6=mr9*sz*tsz;
    f_t tsxg1=tsx*mr1; f_t tsyg4=tsy*mr2; f_t tszg3=tsz*mr3;
    f_t tsxg4=tsx*mr2; f_t tsyg5=tsy*mr5; f_t tszg8=tsz*mr6;
    f_t tsxg7=tsx*mr3; f_t tsyg8=tsy*mr6; f_t tszg9=tsz*mr9;
    f_t rp[3];
    for(unsigned i=0;i<3;i++) {
      rp[i] = static_cast<f_t>(unit_cell.reciprocal_parameters()[i]);
    }
    std::vector<int> mys;
    std::vector<int> mzs;
    f_t cutoff = radius;
    typedef scitbx::math::float_int_conversions<f_t, int> fic;
    for(std::size_t i_site=0;i_site<sites_frac.size();i_site++) {
      fractional<> const& site = sites_frac[i_site];
      f_t xfi=static_cast<f_t>(site[0]);
      f_t yfi=static_cast<f_t>(site[1]);
      f_t zfi=static_cast<f_t>(site[2]);
      f_t cutoffsq=cutoff*cutoff;
      f_t coas = cutoff*rp[0];
      int x1box=fic::ifloor(nx*(xfi-coas));
      int x2box=fic::iceil(nx*(xfi+coas));
      f_t cobs = cutoff*rp[1];
      int y1box=fic::ifloor(ny*(yfi-cobs));
      int y2box=fic::iceil(ny*(yfi+cobs));
      f_t cocs = cutoff*rp[2];
      int z1box=fic::ifloor(nz*(zfi-cocs));
      int z2box=fic::iceil(nz*(zfi+cocs));
      f_t sxbcen=xfi-x1box*sx;
      f_t sybcen=yfi-y1box*sy;
      f_t szbcen=zfi-z1box*sz;
      f_t distsm=mr1*sxbcen*sxbcen+mr5*sybcen*sybcen+mr9*szbcen*szbcen
            +tmr2*sxbcen*sybcen+tmr3*sxbcen*szbcen+tmr6*sybcen*szbcen;
      f_t w7=tsxg1*sxbcen+tsxg4*sybcen+tsxg7*szbcen;
      f_t w8=tsyg4*sxbcen+tsyg5*sybcen+tsyg8*szbcen;
      f_t w9=tszg3*sxbcen+tszg8*sybcen+tszg9*szbcen;
      mys.clear();
      mys.reserve(y2box-y1box+1);
      for (int ky = y1box; ky <= y2box; ky++) {
        int my = ky % ny;
        if (my < 0) my += ny;
        mys.push_back(my);
      }
      mzs.clear();
      mzs.reserve(z2box-z1box+1);
      for (int kz = z1box; kz <= z2box; kz++) {
        int mz = kz % nz;
        if (mz < 0) mz += nz;
        mzs.push_back(mz);
      }
      DataType density_sum = 0;
      std::size_t n_density_contributions = 0;
      f_t distsx = distsm;
      f_t s1xx = sxsq - w7;
      f_t s1xy = sysq - w8;
      f_t s1xz = szsq - w9;
      for (int kx = x1box; kx <= x2box; kx++) {
        int mx = kx % nx;
        if (mx < 0) mx += nx;
        int mxny = mx * ny;
        f_t s2yz = s1xz;
        f_t s2_incr = s1xy;
        f_t s2 = distsx;
        std::vector<int>::const_iterator mye = mys.end();
        for (std::vector<int>::const_iterator
               myi=mys.begin();
               myi!=mye;
               myi++) {
          f_t s3_incr = s2yz;
          f_t dist = s2;
          const DataType* data_mxnymynz = &data[(mxny + (*myi)) * nz];
          std::vector<int>::const_iterator mze = mzs.end();
          for (std::vector<int>::const_iterator
                 mzi=mzs.begin();
                 mzi!=mze;
                 mzi++) {
            if (dist < cutoffsq) {
              density_sum += data_mxnymynz[*mzi];
              n_density_contributions++;
            }
            dist += s3_incr;
            s3_incr += w6;
          }
          s2 += s2_incr;
          s2_incr += w4;
          s2yz += w5;
        }
        distsx += s1xx;
        s1xx += w1;
        s1xy += w2;
        s1xz += w3;
      }
      if (n_density_contributions > 0) {
        result.push_back(  density_sum
                         / static_cast<DataType>(n_density_contributions));
      }
      else {
        result.push_back(0);
      }
    }
    return result;
  }

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_AVERAGE_DENSITIES_H
