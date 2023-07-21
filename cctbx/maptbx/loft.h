#ifndef CCTBX_MAPTBX_LOFT_H
#define CCTBX_MAPTBX_LOFT_H

#include <scitbx/array_family/accessors/c_grid.h>

namespace cctbx { namespace maptbx {

/*

LoFT - Local Fourier Transform.

Authors: Pavel Afonine, Alexandre Urzhumtsev

Release date: mid-July 2023

Highly-optimized local FT. Assumes most map voxels are zeros and only small
subset of voxels are non-zero and so they contribute to summation.

*/

class loft {
public:
  af::shared<std::complex<double> > result;
  af::shared<cctbx::miller::index<> > all_indices;
  int nx, ny, nz;
  int maxh, maxk, maxl;
  af::shared<int> maxkl;
  af::shared<af::shared<int> > maxhkl;
  double two_pi, stepx,stepy,stepz, dcx,dsx, sressq, dtable;
  int maxtab;
  af::shared<double> tabcos;
  af::shared<double> tabsin;
  af::const_ref<int, af::c_grid<3> > const& map_data;

  loft(
    af::const_ref<cctbx::miller::index<> > const& miller_indices,
    af::const_ref<int, af::c_grid<3> > const& map_data_,
    af::shared<double> const& abc,
    double const& d_min)
  :
    nx(map_data_.accessor()[0]),
    ny(map_data_.accessor()[1]),
    nz(map_data_.accessor()[2]),
    two_pi(scitbx::constants::two_pi),
    maxtab(10000),
    sressq(1.0 / (d_min*d_min)),
    map_data(map_data_)
  {
      //
      // Various initializations
      //
      dtable = two_pi / maxtab;
      for(int i = 0; i < maxtab+1; i++) {
        tabcos.push_back(std::cos(i*dtable));
        tabsin.push_back(std::sin(i*dtable));
      }

      maxh = int(abc[0] / d_min) + 1;
      maxk = int(abc[1] / d_min) + 1;
      maxl = int(abc[2] / d_min) + 1;

      maxkl.resize(maxl);
      for (int i = 0; i < maxk; i++) {
        af::shared<int> tmp(maxl);
        maxhkl.push_back(tmp);
      }

      int mhkltot = 0;
      for (int il = 0; il < maxl; il++) {
        double ssl = std::pow(il/abc[2], 2);
        maxkl[il] = int(abc[1] * std::sqrt(sressq - ssl));
        for (int ik = 0; ik < maxkl[il]; ik++) {
           double sskl = std::pow(ik/abc[1], 2) + ssl;
           int tmp = int(abc[0] * std::sqrt(sressq - sskl));
           maxhkl[ik][il] = tmp;
           mhkltot = mhkltot + maxhkl[ik][il];
        }
      }
      result.resize(mhkltot, std::complex<double>(0,0));

      stepx = two_pi / nx;
      stepy = two_pi / ny;
      stepz = two_pi / nz;

      dcx = std::cos(stepx);
      dsx = std::sin(stepx);

      for (int il = 0; il < maxl; il++) {
        for (int ik = 0; ik < maxkl[il]; ik++) {
          for (int ih = 0; ih < maxhkl[ik][il]; ih++) {
            all_indices.push_back(cctbx::miller::index<>(ih,ik,il));
      }}}
      //
      // Main calculations
      //
      double delht  = 0;
      double arg    = 0;
      int karg      = 0;
      double aarg   = 0;
      int iarg      = 0;
      double delarg = 0;
      double delhc  = 0;
      double delkc  = 0;
      double dellc  = 0;
      double delhs  = 0;
      double delks  = 0;
      double dells  = 0;
      double sfrtmp = 0;

      for (int iz = 0; iz < nz; iz++) {
        for (int iy = 0; iy < ny; iy++) {
          int jx = nx;
          for (int ix = 0; ix < nx; ix++) {
            if(map_data(ix,iy,iz)>0) {

              if(ix == jx+1) {
                 delht = delhc * dcx - delhs * dsx;
                 delhs = delhs * dcx + delhc * dsx;
                 delhc = delht;
              }
              else {
                 arg  = ix * stepx;
                 karg = int(arg / two_pi);
                 arg  = arg - karg * two_pi;
                 aarg = arg / dtable;
                 iarg = int(aarg);
                 delarg = aarg - iarg;
                 delhc = tabcos[iarg] * (1.0-delarg) + tabcos[iarg+1] * delarg;
                 delhs = tabsin[iarg] * (1.0-delarg) + tabsin[iarg+1] * delarg;
              }

              jx = ix;

              arg  = iy * stepy;
              karg = int(arg / two_pi);
              arg  = arg - karg * two_pi;
              aarg = arg / dtable;
              iarg = int(aarg);
              delarg = aarg - iarg;
              delkc = tabcos[iarg] * (1.0-delarg) + tabcos[iarg+1] * delarg;
              delks = tabsin[iarg] * (1.0-delarg) + tabsin[iarg+1] * delarg;

              arg   = iz * stepz;
              karg = int(arg / two_pi);
              arg  = arg - karg * two_pi;
              aarg = arg / dtable;
              iarg = int(aarg);
              delarg = aarg - iarg;
              dellc = tabcos[iarg] * (1.0-delarg) + tabcos[iarg+1] * delarg;
              dells = tabsin[iarg] * (1.0-delarg) + tabsin[iarg+1] * delarg;

              //define initial values for the reflection (000)
              //int iref = 0;
              double sfrref = 1.0;
              double sfiref = 0.0;

              double sfrl = sfrref;
              double sfil = sfiref;

              // This is a little bit faster than result[iref]
              std::complex<double>* result_ = result.begin();

              for (int il = 0; il < maxl; il++) {
                 double sfrkl = sfrl;
                 double sfikl = sfil;
                 for (int ik = 0; ik < maxkl[il]; ik++) {
                    double sfrhkl = sfrkl;
                    double sfihkl = sfikl;
                    for (int ih = 0; ih < maxhkl[ik][il]; ih++) {
                       //result[iref] += std::complex<double>(sfrhkl,sfihkl);
                       *result_ += std::complex<double>(sfrhkl,sfihkl);
                       ++result_;

                       //iref = iref + 1;
                       sfrtmp = sfrhkl * delhc - sfihkl * delhs;
                       sfihkl = sfihkl * delhc + sfrhkl * delhs;
                       sfrhkl = sfrtmp;
                    }
                    sfrtmp = sfrkl * delkc - sfikl * delks;
                    sfikl  = sfikl * delkc + sfrkl * delks;
                    sfrkl  = sfrtmp;
                 }
                 sfrtmp = sfrl * dellc - sfil * dells;
                 sfil   = sfil * dellc + sfrl * dells;
                 sfrl   = sfrtmp;
              }
            }
      }}}
    }

  af::shared<std::complex<double> >  structure_factors() { return result; }
  af::shared<cctbx::miller::index<> > indices() { return all_indices; }

};


// Straightforward reference implementation using the formula directly

//af::shared<std::complex<double> > sf_from_map_simple(
//  af::const_ref<cctbx::miller::index<> > const& miller_indices,
//  af::const_ref<int, af::c_grid<3> > const& map_data,
//  af::const_ref<double> const& sin_table_,
//  af::const_ref<double> const& cos_table_,
//  double const& step,
//  int const& n)
//{
//
//    math::cos_sin_table<double> cos_sin(1000);
//
//    int hkl_size = miller_indices.size();
//    af::shared<std::complex<double> > result(hkl_size);
//    af::shared<scitbx::vec3<int> > jxyz;
//    int nx = map_data.accessor()[0];
//    int ny = map_data.accessor()[1];
//    int nz = map_data.accessor()[2];
//    for (int jx = 0; jx < nx; jx++) {
//      for (int jy = 0; jy < ny; jy++) {
//        for (int jz = 0; jz < nz; jz++) {
//          if(map_data(jx,jy,jz)==1) jxyz.push_back(scitbx::vec3<int>(jx,jy,jz));
//    }}}
//    af::shared<double> hi(hkl_size);
//    af::shared<double> ki(hkl_size);
//    af::shared<double> li(hkl_size);
//    double pi = scitbx::constants::pi;
//    double tponx = 2*pi/nx, tpony = 2*pi/ny, tponz = 2*pi/nz;
//    for (int i = 0; i < hkl_size; i++) {
//      cctbx::miller::index<> mi = miller_indices[i];
//      hi[i] = tponx*mi[0]/scitbx::constants::two_pi;
//      ki[i] = tpony*mi[1]/scitbx::constants::two_pi;
//      li[i] = tponz*mi[2]/scitbx::constants::two_pi;
//    }
//    int reg_size = jxyz.size();
//    for (int i=0; i<hkl_size; i++) {
//      double h = hi[i];
//      double k = ki[i];
//      double l = li[i];
//      std::complex<double> arg(0,0);
//      for (int j=0; j<reg_size; j++) {
//        scitbx::vec3<int> p = jxyz[j];
//        double arg_ = h*p[0] + k*p[1] + l*p[2];
//        // Use standard sin and cos
//        //arg += std::complex<double>(std::cos(arg_), std::sin(arg_));
//        // Use sin/cos tables
//        //double co = scitbx::math::cos_table(cos_table_,arg_,step,n,false);
//        //double si = scitbx::math::sin_table(sin_table_,arg_,step,n,false);
//        //arg += std::complex<double>(co, si);
//        arg += cos_sin.get(arg_);
//
//      }
//      result[i]=arg;
//    }
//    return result;
//}

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_LOFT_H
