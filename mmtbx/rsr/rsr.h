#ifndef MMTBX_RSR_H
#define MMTBX_RSR_H

#include <cctbx/sgtbx/space_group.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <vector>
#include <scitbx/array_family/shared.h>
#include <cctbx/import_scitbx_af.h>
#include <cmath>
#include <cctbx/adptbx.h>
#include <assert.h>
#include <math.h>
#include <iostream>
#include <mmtbx/error.h>
#include <cctbx/xray/scattering_type_registry.h>

#include <cctbx/xray/sampling_base.h>

#include <cctbx/maptbx/bcr/bcr.h>
#include <boost/python/list.hpp>

using namespace std;
namespace mmtbx { namespace rsr {
namespace af=scitbx::af;
using scitbx::vec3;
using scitbx::mat3;
using scitbx::sym_mat3;

template <typename FloatType=double,
          typename XrayScattererType=cctbx::xray::scatterer<> >
class gaussian_density {
public:
  typedef FloatType f_t;

  gaussian_density(
     XrayScattererType scatterer,
     cctbx::xray::scattering_type_registry const& scattering_type_registry,
     cctbx::xray::detail::exponent_table<f_t>& exp_table,
     f_t wing_cutoff)
  :
  exp_table_(&exp_table), wing_cutoff_(wing_cutoff)
  {
   cctbx::eltbx::xray_scattering::gaussian const& gaussian =
     scattering_type_registry.gaussian_not_optional(scatterer.scattering_type);
   f_t b_iso = cctbx::adptbx::u_as_b(scatterer.u_iso);
   n_terms = gaussian.n_terms();
   for(std::size_t i=0;i<n_terms;i++) {
     scitbx::math::gaussian::term<f_t> gti = gaussian.terms()[i];
     f_t b_all = gti.b + b_iso;
     f_t d = b_all * b_all * b_all;
     MMTBX_ASSERT(d > 0.);
     as[i] = cctbx::xray::detail::eight_pi_pow_3_2 * scatterer.weight() *
             gti.a / std::sqrt(d);
     bs[i] = -cctbx::xray::detail::four_pi_sq / b_all;
   }
  }

  f_t atom_radius(f_t max_radius=5.0, f_t step=0.1) const
  {
    f_t rho_0(0);
    for(std::size_t i=0;i<n_terms;i++) {
      rho_0 += as[i];
    }
    f_t radius(0);
    f_t rho(0);
    while(radius < max_radius) {
      f_t rho(0);
      f_t r_sq = radius * radius;
      for(std::size_t i=0;i<n_terms;i++) {
        rho += as[i]*exp_table(bs[i] * r_sq);
      }
      if(rho/rho_0 < wing_cutoff_) return radius;
      radius += step;
    }
    return radius;
  }

  f_t rho(f_t const& d_sq) const
  {
    f_t r(0);
    for(std::size_t i=0;i<n_terms;i++) {
      //r += as[i]*std::exp(bs[i]*d_sq); // XXX slow but exact
      r += as[i]*exp_table(bs[i] * d_sq);
    }
    return r;
  }

  f_t exp_table(f_t const& x) const { return (*exp_table_)(x); }

protected:
  BOOST_STATIC_CONSTANT(std::size_t, max_rho_terms =
    cctbx::eltbx::xray_scattering::gaussian::max_n_terms+1);
  af::tiny<f_t, max_rho_terms> as;
  af::tiny<f_t, max_rho_terms> bs;
  std::size_t n_terms;
  cctbx::xray::detail::exponent_table<f_t>* exp_table_;
  f_t wing_cutoff_;
};


template <typename FloatType=double,
          typename XrayScattererType=cctbx::xray::scatterer<> >
class manager {
public:
  typedef FloatType f_t;
  af::versa<f_t, af::c_grid<3> > density_array;
  manager(int const& nx,
          int const& ny,
          int const& nz,
          cctbx::xray::scattering_type_registry const& scattering_type_registry,
          cctbx::uctbx::unit_cell const& unit_cell,
          af::const_ref<XrayScattererType> const& scatterers,
          f_t const& exp_table_one_over_step_size=-100,
          f_t const& wing_cutoff=1.e-3)
  {
    f_t dist=0.0, dist_=0.0;
    //af::versa<f_t, af::c_grid<3> > density_array(af::c_grid<3>(nx,ny,nz),
    //    af::init_functor_null<f_t>());
    density_array.resize(af::c_grid<3>(nx,ny,nz), 0);
    //
    scitbx::sym_mat3<f_t> metrical_matrix = unit_cell.metrical_matrix();
    f_t mr1= static_cast<f_t>(metrical_matrix[0]);
    f_t mr5= static_cast<f_t>(metrical_matrix[1]);
    f_t mr9= static_cast<f_t>(metrical_matrix[2]);
    f_t mr2= static_cast<f_t>(metrical_matrix[3]);
    f_t mr3= static_cast<f_t>(metrical_matrix[4]);
    f_t mr6= static_cast<f_t>(metrical_matrix[5]);
    f_t tmr2=mr2*2; f_t sx=1/static_cast<f_t>(nx); f_t tsx=sx*2;
    f_t tmr3=mr3*2; f_t sy=1/static_cast<f_t>(ny); f_t tsy=sy*2;
    f_t tmr6=mr6*2; f_t sz=1/static_cast<f_t>(nz); f_t tsz=sz*2;
    f_t sxsq=mr1*sx*sx; f_t w1=mr1*sx*tsx; f_t w4=mr5*sy*tsy;
    f_t sysq=mr5*sy*sy; f_t w2=mr2*sx*tsy; f_t w5=mr6*sy*tsz;
    f_t szsq=mr9*sz*sz; f_t w3=mr3*sx*tsz; f_t w6=mr9*sz*tsz;
    f_t tsxg1=tsx*mr1; f_t tsyg4=tsy*mr2; f_t tszg3=tsz*mr3;
    f_t tsxg4=tsx*mr2; f_t tsyg5=tsy*mr5; f_t tszg8=tsz*mr6;
    f_t tsxg7=tsx*mr3; f_t tsyg8=tsy*mr6; f_t tszg9=tsz*mr9;
    f_t rp[3];
    for(unsigned i=0;i<3;i++) {
      rp[i] = static_cast<f_t>(unit_cell.reciprocal_parameters()[i]);
    }
    //
    cctbx::xray::detail::exponent_table<f_t>
      exp_table(exp_table_one_over_step_size);
    //
    f_t* density_array_ = density_array.begin();
    for(std::size_t i_site=0;i_site<scatterers.size();i_site++) {
      XrayScattererType const& scatterer = scatterers[i_site];
      gaussian_density<FloatType, XrayScattererType> gd(scatterer,
        scattering_type_registry, exp_table, wing_cutoff);
      f_t cutoff=gd.atom_radius();
      f_t cutoffsq=cutoff*cutoff;
      f_t coas = cutoff*rp[0];
      f_t cobs = cutoff*rp[1];
      f_t cocs = cutoff*rp[2];
      cctbx::fractional<> const& site = scatterer.site;
      f_t xfi= static_cast<f_t>(site[0]);
      f_t yfi= static_cast<f_t>(site[1]);
      f_t zfi= static_cast<f_t>(site[2]);
      int x1box=ifloor(nx*(xfi-coas));
      int x2box= iceil(nx*(xfi+coas));
      int y1box=ifloor(ny*(yfi-cobs));
      int y2box= iceil(ny*(yfi+cobs));
      int z1box=ifloor(nz*(zfi-cocs));
      int z2box= iceil(nz*(zfi+cocs));
      f_t sxbcen=xfi-x1box*sx;
      f_t sybcen=yfi-y1box*sy;
      f_t szbcen=zfi-z1box*sz;
      f_t distsm=mr1*sxbcen*sxbcen+mr5*sybcen*sybcen+mr9*szbcen*szbcen
            +tmr2*sxbcen*sybcen+tmr3*sxbcen*szbcen+tmr6*sybcen*szbcen;
      f_t w7=tsxg1*sxbcen+tsxg4*sybcen+tsxg7*szbcen;
      f_t w8=tsyg4*sxbcen+tsyg5*sybcen+tsyg8*szbcen;
      f_t w9=tszg3*sxbcen+tszg8*sybcen+tszg9*szbcen;
      f_t distsx = distsm;
      f_t s1xx = sxsq - w7;
      f_t s1xy = sysq - w8;
      f_t s1xz = szsq - w9;
      for(int kx = x1box; kx <= x2box; kx++) {
        int mx = scitbx::math::mod_positive(kx,nx);
        int mxny = mx*ny;
        //f_t xn=xfi-double(kx)/nx; //XXX
        f_t s2yz = s1xz;
        f_t s2_incr = s1xy;
        f_t s2 = distsx;
        for(int ky = y1box; ky <= y2box; ky++) {
          int my = scitbx::math::mod_positive(ky,ny);
          int mxnypmynz = (mxny+my)*nz;
          //f_t yn=yfi-double(ky)/ny; //XXX
          f_t s3_incr = s2yz;
          f_t dist = s2;
          for (int kz = z1box; kz <= z2box; kz++) {
            //f_t zn=zfi-double(kz)/nz; //XXX
            //dist_=mr1*xn*xn+mr5*yn*yn+mr9*zn*zn+tmr2*xn*yn+tmr3*xn*zn+tmr6*yn*zn;
            //MMTBX_ASSERT(std::abs(dist - dist_) < 1.e-3);
            if(dist <= cutoffsq) {
              int mz = scitbx::math::mod_positive(kz,nz);
              density_array_[mxnypmynz+mz] += gd.rho(dist); // (mx*NY+my)*NZ+mz
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
    }
  }

protected:
  inline static int ifloor(f_t const& x)
  {
    return scitbx::math::float_int_conversions<f_t, int>::ifloor(x);
  }
  inline static int iceil(f_t const& x)
  {
    return scitbx::math::float_int_conversions<f_t, int>::iceil(x);
  }
};

//----------------------BCR manager --------------------------------------------

template <typename FloatType=double,
          typename XrayScattererType=cctbx::xray::scatterer<>,
          typename BCRSType=cctbx::maptbx::bcr_model<> >
class manager_BCR {
public:
  typedef FloatType f_t;
  af::versa<f_t, af::c_grid<3> > density_array;
  manager_BCR(int const& nx,
          int const& ny,
          int const& nz,
          cctbx::xray::scattering_type_registry const& scattering_type_registry,
          cctbx::uctbx::unit_cell const& unit_cell,
          //af::const_ref<XrayScattererType> const& scatterers,
          boost::python::list const& BCRscatterers,

          f_t const& exp_table_one_over_step_size=-100,
          f_t const& wing_cutoff=1.e-3)
  {
    f_t dist=0.0, dist_=0.0;
    //af::versa<f_t, af::c_grid<3> > density_array(af::c_grid<3>(nx,ny,nz),
    //    af::init_functor_null<f_t>());
    density_array.resize(af::c_grid<3>(nx,ny,nz), 0);
    //
    scitbx::sym_mat3<f_t> metrical_matrix = unit_cell.metrical_matrix();
    f_t mr1= static_cast<f_t>(metrical_matrix[0]);
    f_t mr5= static_cast<f_t>(metrical_matrix[1]);
    f_t mr9= static_cast<f_t>(metrical_matrix[2]);
    f_t mr2= static_cast<f_t>(metrical_matrix[3]);
    f_t mr3= static_cast<f_t>(metrical_matrix[4]);
    f_t mr6= static_cast<f_t>(metrical_matrix[5]);
    f_t tmr2=mr2*2; f_t sx=1/static_cast<f_t>(nx); f_t tsx=sx*2;
    f_t tmr3=mr3*2; f_t sy=1/static_cast<f_t>(ny); f_t tsy=sy*2;
    f_t tmr6=mr6*2; f_t sz=1/static_cast<f_t>(nz); f_t tsz=sz*2;
    f_t sxsq=mr1*sx*sx; f_t w1=mr1*sx*tsx; f_t w4=mr5*sy*tsy;
    f_t sysq=mr5*sy*sy; f_t w2=mr2*sx*tsy; f_t w5=mr6*sy*tsz;
    f_t szsq=mr9*sz*sz; f_t w3=mr3*sx*tsz; f_t w6=mr9*sz*tsz;
    f_t tsxg1=tsx*mr1; f_t tsyg4=tsy*mr2; f_t tszg3=tsz*mr3;
    f_t tsxg4=tsx*mr2; f_t tsyg5=tsy*mr5; f_t tszg8=tsz*mr6;
    f_t tsxg7=tsx*mr3; f_t tsyg8=tsy*mr6; f_t tszg9=tsz*mr9;
    f_t rp[3];
    for(unsigned i=0;i<3;i++) {
      rp[i] = static_cast<f_t>(unit_cell.reciprocal_parameters()[i]);
    }
    //
    cctbx::xray::detail::exponent_table<f_t>
      exp_table(exp_table_one_over_step_size);
    //
    f_t* density_array_ = density_array.begin();
    //for(std::size_t i_site=0;i_site<scatterers.size();i_site++) {
    for(std::size_t i_site=0;i_site<boost::python::len(BCRscatterers);i_site++) {
      //XrayScattererType const& scatterer = scatterers[i_site];
      //BCRSType const& BCRscatterer = BCRscatterers[i_site];

      BCRSType const& BCRscatterer =
         boost::python::extract<BCRSType>(BCRscatterers[i_site])();

      //gaussian_density<FloatType, XrayScattererType> gd(scatterer,
      //  scattering_type_registry, exp_table, wing_cutoff);

      cctbx::maptbx::calculator<f_t> calc =
        cctbx::maptbx::calculator<f_t>(BCRscatterer, exp_table);

      //cctbx::maptbx::calculator<f_t> calc =
      //  cctbx::maptbx::calculator<f_t>(BCRscatterer);

      f_t cutoff=calc.atom_radius();
      f_t cutoffsq=cutoff*cutoff;
      f_t coas = cutoff*rp[0];
      f_t cobs = cutoff*rp[1];
      f_t cocs = cutoff*rp[2];
      cctbx::fractional<> const& site = BCRscatterer.scatterer.site; //
      f_t xfi= static_cast<f_t>(site[0]);
      f_t yfi= static_cast<f_t>(site[1]);
      f_t zfi= static_cast<f_t>(site[2]);
      int x1box=ifloor(nx*(xfi-coas));
      int x2box= iceil(nx*(xfi+coas));
      int y1box=ifloor(ny*(yfi-cobs));
      int y2box= iceil(ny*(yfi+cobs));
      int z1box=ifloor(nz*(zfi-cocs));
      int z2box= iceil(nz*(zfi+cocs));
      f_t sxbcen=xfi-x1box*sx;
      f_t sybcen=yfi-y1box*sy;
      f_t szbcen=zfi-z1box*sz;
      f_t distsm=mr1*sxbcen*sxbcen+mr5*sybcen*sybcen+mr9*szbcen*szbcen
            +tmr2*sxbcen*sybcen+tmr3*sxbcen*szbcen+tmr6*sybcen*szbcen;
      f_t w7=tsxg1*sxbcen+tsxg4*sybcen+tsxg7*szbcen;
      f_t w8=tsyg4*sxbcen+tsyg5*sybcen+tsyg8*szbcen;
      f_t w9=tszg3*sxbcen+tszg8*sybcen+tszg9*szbcen;
      f_t distsx = distsm;
      f_t s1xx = sxsq - w7;
      f_t s1xy = sysq - w8;
      f_t s1xz = szsq - w9;
      for(int kx = x1box; kx <= x2box; kx++) {
        int mx = scitbx::math::mod_positive(kx,nx);
        int mxny = mx*ny;
        //f_t xn=xfi-double(kx)/nx; //XXX
        f_t s2yz = s1xz;
        f_t s2_incr = s1xy;
        f_t s2 = distsx;
        for(int ky = y1box; ky <= y2box; ky++) {
          int my = scitbx::math::mod_positive(ky,ny);
          int mxnypmynz = (mxny+my)*nz;
          //f_t yn=yfi-double(ky)/ny; //XXX
          f_t s3_incr = s2yz;
          f_t dist = s2;
          for (int kz = z1box; kz <= z2box; kz++) {
            //f_t zn=zfi-double(kz)/nz; //XXX
            //dist_=mr1*xn*xn+mr5*yn*yn+mr9*zn*zn+tmr2*xn*yn+tmr3*xn*zn+tmr6*yn*zn;
            //MMTBX_ASSERT(std::abs(dist - dist_) < 1.e-3);
            if(dist <= cutoffsq) {
              int mz = scitbx::math::mod_positive(kz,nz);
              //density_array_[mxnypmynz+mz] += calc.rho(dist); // (mx*NY+my)*NZ+mz

              if(std::abs(dist)<1.e-6) dist=0; // cause of: exponent_table: excessive range. crash

              density_array_[mxnypmynz+mz] += calc.rho(std::sqrt(dist));
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
    }
  }

protected:
  inline static int ifloor(f_t const& x)
  {
    return scitbx::math::float_int_conversions<f_t, int>::ifloor(x);
  }
  inline static int iceil(f_t const& x)
  {
    return scitbx::math::float_int_conversions<f_t, int>::iceil(x);
  }
};

//----------------------BCR manager END ----------------------------------------


}} // namespace mmtbx::rsr

#endif // MMTBX_RSR_H
