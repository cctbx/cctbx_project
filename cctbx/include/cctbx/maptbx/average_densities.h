// Based on code by Pavel Afonine
#ifndef CCTBX_MAPTBX_AVERAGE_DENSITIES_H
#define CCTBX_MAPTBX_AVERAGE_DENSITIES_H

#include <cctbx/uctbx.h>
#include <scitbx/array_family/accessors/c_grid.h>

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
