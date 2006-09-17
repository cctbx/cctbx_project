// Based on code by Pavel Afonine
#ifndef CCTBX_MAPTBX_AVERAGE_DENSITIES_H
#define CCTBX_MAPTBX_AVERAGE_DENSITIES_H

#include <cctbx/uctbx.h>
#include <scitbx/array_family/accessors/c_grid.h>

using scitbx::vec3;
namespace af=scitbx::af;

namespace cctbx { namespace maptbx {

//! Fortran-like nearest integer.
inline
int
nint(double x)
{
  return int(std::ceil(x+0.5)-(std::fmod(x*0.5+0.25,1.0)!=0));
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
                    mx = cctbx::math::mod_positive(kx, NX);
                    my = cctbx::math::mod_positive(ky, NY);
                    mz = cctbx::math::mod_positive(kz, NZ);
                    r_ave_xyz += data(mx,my,mz);
                    counter += 1;
            }}}
            new_data_ref(lx,ly,lz) = r_ave_xyz / counter;
    }}}
    return new_data;
}


class grid_points_in_sphere_around_atom_and_distances {
public:
  grid_points_in_sphere_around_atom_and_distances(
                           cctbx::uctbx::unit_cell const& uc,
                           af::const_ref<double, af::c_grid<3> > const& data,
                           double const& rad,
                           double const& shell,
                           vec3<double> const& site_frac)
  {
    int nx = data.accessor()[0];
    int ny = data.accessor()[1];
    int nz = data.accessor()[2];
    double distsq;
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
       int x1box = nint( nx*(xf-coas) ) - 1;
       int x2box = nint( nx*(xf+coas) ) + 1;
       int y1box = nint( ny*(yf-cobs) ) - 1;
       int y2box = nint( ny*(yf+cobs) ) + 1;
       int z1box = nint( nz*(zf-cocs) ) - 1;
       int z2box = nint( nz*(zf+cocs) ) + 1;
       for(kx = x1box; kx <= x2box; kx++) {
       //for(kx = 0; kx < nx; kx++) {
           double xn=xf-double(kx)/nx;
           for(ky = y1box; ky <= y2box; ky++) {
           //for(ky = 0; ky < ny; ky++) {
               double yn=yf-double(ky)/ny;
               for(kz = z1box; kz <= z2box; kz++) {
               //for(kz = 0; kz < nz; kz++) {
                   double zn=zf-double(kz)/nz;
                   distsq=mr1*xn*xn+mr5*yn*yn+mr9*zn*zn+
                                              tmr2*xn*yn+tmr3*xn*zn+tmr6*yn*zn;
                   if(distsq <= radsq) {
                   //if(1) {
                      mx = cctbx::math::mod_positive(kx, nx);
                      my = cctbx::math::mod_positive(ky, ny);
                      mz = cctbx::math::mod_positive(kz, nz);
                      //double rho = 2.0 * std::exp (- 3.0 * distsq );
                      data_at_grid_points_.push_back(data(mx,my,mz));
                      //data_at_grid_points_.push_back(rho);
                      distances_.push_back(std::sqrt(distsq));
                   }
               }
           }
       }
    }
  }
  af::shared<double> data_at_grid_points() { return data_at_grid_points_; }
  af::shared<double> distances() { return distances_; }
protected:
  af::shared<double> data_at_grid_points_;
  af::shared<double> distances_;
};

class one_gaussian_peak_approximation {
public:
  one_gaussian_peak_approximation(
                              af::const_ref<double> const& data_at_grid_points,
                              af::const_ref<double> const& distances)
  {
           double p=0.0, q=0.0, r=0.0, s=0.0;
           int n = 0;
           int number_of_points = data_at_grid_points.size();
           for(int i = 0; i < number_of_points; i++) {
            if(distances[i] > 1.e-4 && data_at_grid_points[i] > 0.0) {
             n = n + 1;
              double data_i = data_at_grid_points[i];
              double distance_i = distances[i];
              double distance_i_sq = distance_i*distance_i;
               p=p+std::log(data_i);
               q=q+distance_i_sq;
               r=r+(distance_i_sq*distance_i_sq);
               s=s+distance_i_sq*std::log(data_i);
      }
    }
    double alpha_opt = (p-s*q/r) / (double(n)-q*q/r);
    a_real_space_ = std::exp( alpha_opt );
           b_real_space_ = 1/r*(alpha_opt*q - s);

           double tmp = b_real_space_/scitbx::constants::pi;
           double pi_sq = scitbx::constants::pi*scitbx::constants::pi;
           a_reciprocal_space_ = a_real_space_/std::sqrt(tmp*tmp*tmp);
           b_reciprocal_space_ = pi_sq/b_real_space_*4; //      not deja divise par 4 !!!
           double num = 0.0;
           double denum = 0.0;
           for(int i = 0; i < number_of_points; i++) {
            if(distances[i] > 1.e-4 && data_at_grid_points[i] > 0.0) {
              double data_i = data_at_grid_points[i];
              double distance_i = distances[i];
              double distance_i_sq = distance_i*distance_i;
              double data_i_approx =
                              a_real_space_ * std::exp(-b_real_space_*distance_i_sq);
               num=num+std::abs(data_i-data_i_approx);
               denum=denum+data_i_approx;
      }
    }
    gof_ = num/denum*100.;
  }
  double a_real_space() { return a_real_space_; }
  double b_real_space() { return b_real_space_; }
  double a_reciprocal_space() { return a_reciprocal_space_; }
  double b_reciprocal_space() { return b_reciprocal_space_; }
  double gof() { return gof_; }
protected:
  double a_real_space_;
  double b_real_space_;
  double a_reciprocal_space_;
  double b_reciprocal_space_;
  double gof_;
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
