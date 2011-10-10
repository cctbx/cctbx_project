// Based on code by Pavel Afonine
#ifndef CCTBX_MAPTBX_AVERAGE_DENSITIES_H
#define CCTBX_MAPTBX_AVERAGE_DENSITIES_H

#include <cctbx/uctbx.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/math/utils.h>
#include <cctbx/maptbx/eight_point_interpolation.h>

namespace cctbx { namespace maptbx {

//! Fortran-like nearest integer.
inline
int
nint(double x)
{
  return int(std::ceil(x+0.5)-(std::fmod(x*0.5+0.25,1.0)!=0));
}

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
                   scitbx::vec3<double> const& translation_vector)
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
      for (int k = 0; k < nz; k++) {
        cctbx::fractional<> grid_frac_in_1 = cctbx::fractional<>(gpfx, gpfy,
          k/static_cast<double>(nz));
        cctbx::cartesian<> grid_cart_in_1 = unit_cell_2.orthogonalize(
          grid_frac_in_1);
        cctbx::cartesian<> point_cart_in_2 = rotation_matrix * grid_cart_in_1 +
          translation_vector;
        cctbx::fractional<> point_frac_in_2 =
          unit_cell_1.fractionalize(point_cart_in_2);
        result_map_ref(i,j,k) = tricubic_interpolation(map_data_1,
          point_frac_in_2);
  }}}
  return result_map;
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

af::versa<double, af::c_grid<3> > box_map_averaging(
                             cctbx::uctbx::unit_cell const& unit_cell,
                             af::const_ref<double, af::c_grid<3> > const& data,
                             double const& rad)
{
    int lx,ly,lz,kx,ky,kz,mx,my,mz,x1box,x2box,y1box,y2box,z1box,z2box;
    int NX = data.accessor()[0];
    int NY = data.accessor()[1];
    int NZ = data.accessor()[2];
    double xrad = rad*unit_cell.reciprocal_parameters()[0]*NX;
    double yrad = rad*unit_cell.reciprocal_parameters()[1]*NY;
    double zrad = rad*unit_cell.reciprocal_parameters()[2]*NZ;
    af::versa<double, af::c_grid<3> > new_data;
    new_data.resize(af::c_grid<3>(NX,NY,NZ), 0.0);
    af::ref<double, af::c_grid<3> > new_data_ref = new_data.ref();

    for (lx = 0; lx < NX; lx++) {
      for (ly = 0; ly < NY; ly++) {
        for (lz = 0; lz < NZ; lz++) {
            double r_ave_xyz = 0.0;
            int counter = 0;
            x1box=nint(static_cast<double>(lx)-xrad) - 1;
            x2box=nint(static_cast<double>(lx)+xrad) + 1;
            y1box=nint(static_cast<double>(ly)-yrad) - 1;
            y2box=nint(static_cast<double>(ly)+yrad) + 1;
            z1box=nint(static_cast<double>(lz)-zrad) - 1;
            z2box=nint(static_cast<double>(lz)+zrad) + 1;
            for (kx = x1box; kx <= x2box; kx++) {
              for (ky = y1box; ky <= y2box; ky++) {
                for (kz = z1box; kz <= z2box; kz++) {
                    mx = scitbx::math::mod_positive(kx, NX);
                    my = scitbx::math::mod_positive(ky, NY);
                    mz = scitbx::math::mod_positive(kz, NZ);
                    r_ave_xyz += data(mx,my,mz);
                    counter += 1;
            }}}
            new_data_ref(lx,ly,lz) = r_ave_xyz / counter;
    }}}
    return new_data;
}

class cumulative_histogramm {
public:
  af::shared<double> sigs;
  af::shared<double> h1;
  af::shared<double> h2;
  af::shared<double> h12;
  cumulative_histogramm(
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
    double max1 = af::max(map_1);
    double min1 = af::min(map_1);
    double max2 = af::max(map_2);
    double min2 = af::min(map_2);
    double start = std::min(min1,min2);
    double end   = std::max(max1,max2);
    //double hsize = 100./map_1.size();
    double hsize = 1./map_1.size();
    double inc = (end-start)/100.;
    double sigp = 0;
    double sigm = 0;
    for (double sig=start; sig<=end; sig += inc) {
      af::tiny<int, 2> r = map_loop(map_1, map_2, sig);
      double h1_ = r[0]*hsize;
      double h2_ = r[1]*hsize;
      if(h1_<0.01 && h2_<0.01 && sigp==0)  sigp = sig;
      if(h1_>0.99 && h2_>0.99) sigm = sig;
    }
    inc = (sigm-start)/100;
    for (double sig=start; sig<sigm; sig += inc) {
      af::tiny<int, 2> r = map_loop(map_1, map_2, sig);
      sigs.push_back(sig);
      double h1_ = r[0]*hsize;
      double h2_ = r[1]*hsize;
      h1.push_back(h1_);
      h2.push_back(h2_);
      h12.push_back((h1_+h2_)/2);
    }
    inc = (sigp-sigm)/100.;
    for (double sig=sigm; sig<sigp; sig += inc) {
      af::tiny<int, 2> r = map_loop(map_1, map_2, sig);
      sigs.push_back(sig);
      double h1_ = r[0]*hsize;
      double h2_ = r[1]*hsize;
      h1.push_back(h1_);
      h2.push_back(h2_);
      h12.push_back((h1_+h2_)/2);
    }
    inc = (end-sigp)/100.;
    for (double sig=sigp; sig<=end+1.e-6*end; sig += inc) {
      af::tiny<int, 2> r = map_loop(map_1, map_2, sig);
      sigs.push_back(sig);
      double h1_ = r[0]*hsize;
      double h2_ = r[1]*hsize;
      h1.push_back(h1_);
      h2.push_back(h2_);
      h12.push_back((h1_+h2_)/2);
    }
  }

  af::shared<double> histogram_1() {return h1;}
  af::shared<double> histogram_2() {return h2;}
  af::shared<double> histogram_average() {return h12;}
  af::shared<double> values() {return sigs;}

protected:
  af::tiny<int, 2> map_loop(
    af::const_ref<double, af::c_grid<3> > m1,
    af::const_ref<double, af::c_grid<3> > m2,
    double sig)
  {
    int nx = m1.accessor()[0];
    int ny = m1.accessor()[1];
    int nz = m1.accessor()[2];
    int p1=0, p2=0;
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
          if(m1(i,j,k) > sig) p1 += 1;
          if(m2(i,j,k) > sig) p2 += 1;
        }
      }
    }
    return af::tiny<int, 2> (p1,p2);
  }
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
          //map_2_new(i,j,k)=map_2(i,j,k);
          //map_2_new(i,j,k)=(map_2(i,j,k)-min2)*(max1-min1)/(max2-min2)+min1;
          //map_1_new(i,j,k)=map_1(i,j,k);
          double m2 = (map_2(i,j,k)-min2)*(2)/(max2-min2)-1;
          double m1 = (map_1(i,j,k)-min1)*(2)/(max1-min1)-1;
          map_2_new(i,j,k)=m2;
          map_1_new(i,j,k)=m1;
          mean1 = mean1 + m1;
          mean2 = mean2 + m2;
        }
      }
    }
    CCTBX_ASSERT(std::abs(af::max(map_1_new.ref())-af::max(map_2_new.ref()))<1e-6);
    CCTBX_ASSERT(std::abs(af::min(map_1_new.ref())-af::min(map_2_new.ref()))<1e-6);
    mean1 = mean1/map_1.size();
    mean2 = mean2/map_1.size();
    for (int i = 0; i < nx1; i++) {
      for (int j = 0; j < ny1; j++) {
        for (int k = 0; k < nz1; k++) {
          map_2_new(i,j,k)=map_2_new(i,j,k) - mean2;
          map_1_new(i,j,k)=map_1_new(i,j,k) - mean1;
        }
      }
    }
    cumulative_histogramm ch = cumulative_histogramm(map_1_new.ref(),map_2_new.ref());
    h1 = ch.histogram_1();
    h2 = ch.histogram_2();
    h12 = ch.histogram_average();
    sigs = ch.values();
    //
    map_mod(map_1_new.ref(), h1.ref(), sigs.ref(), h12.ref());
    map_mod(map_2_new.ref(), h2.ref(), sigs.ref(), h12.ref());
  }

  af::versa<double, af::c_grid<3> > map_1() {return map_1_new;}
  af::versa<double, af::c_grid<3> > map_2() {return map_2_new;}

protected:
  void map_mod(
    af::ref<double, af::c_grid<3> > m,
    af::const_ref<double> h,
    af::const_ref<double> sigs,
    af::const_ref<double> h_ave) {
    int nx = m.accessor()[0];
    int ny = m.accessor()[1];
    int nz = m.accessor()[2];
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
          //
          double m2n = m(i,j,k);
          double sig2_best = sigs[0];
          double r = 1e+9;
          double h2_best = 0;
          for (int n=0; n<sigs.size(); n++) {
            double r_ = std::abs(m2n-sigs[n]);
            if(r_<r) {
              r = r_;
              sig2_best = sigs[n];
              h2_best = h[n];
            }
          }
          //
          double sig1_best = sigs[0];
          double h1_best = 0;
          r = 1e+9;
          for (int n=0; n<h1.size(); n++) {
            double r_ = std::abs(h2_best-h_ave[n]);
            if(r_<r) {
              r = r_;
              sig1_best = sigs[n];
              h1_best = h_ave[n];
            }
          }
          m(i,j,k)=sig1_best;
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
