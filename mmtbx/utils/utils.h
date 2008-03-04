#ifndef MMTBX_UTILS_H
#define MMTBX_UTILS_H

#include <scitbx/array_family/shared.h>
#include <mmtbx/error.h>
#include <cctbx/uctbx.h>
//
#include <cctbx/xray/scatterer.h>
#include <cctbx/sgtbx/sym_equiv_sites.h>
#include <cctbx/sgtbx/site_symmetry_table.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <cstdio>
//

using namespace std;
namespace mmtbx { namespace utils {
namespace af=scitbx::af;
using scitbx::mat3;
using cctbx::uctbx::unit_cell;

template <typename FloatType=double, typename cctbx_frac=cctbx::fractional<> >
class fit_hoh
{
  public:
    cctbx_frac site_cart_o_fitted;
    cctbx_frac site_cart_h1_fitted;
    cctbx_frac site_cart_h2_fitted;
    scitbx::vec3<FloatType> origin_cart;
    FloatType dist_best;

    fit_hoh() {}

    fit_hoh(cctbx_frac const& site_frac_o,
            cctbx_frac const& site_frac_h1,
            cctbx_frac const& site_frac_h2,
            cctbx_frac const& site_frac_peak1,
            cctbx_frac const& site_frac_peak2,
            FloatType const& angular_shift,
            cctbx::uctbx::unit_cell const& unit_cell)
    :
    origin_cart(unit_cell.orthogonalize(site_frac_o)),
    site_cart_o_fitted(unit_cell.orthogonalize(site_frac_o)),
    site_cart_h1_fitted(unit_cell.orthogonalize(site_frac_h1)),
    site_cart_h2_fitted(unit_cell.orthogonalize(site_frac_h2)),
    dist_best(1.e+9)
    {
      CCTBX_ASSERT(angular_shift > 0 && angular_shift < 360);
      bool is_one_peak = false;
      FloatType diff = unit_cell.distance_sq(site_frac_peak1, site_frac_peak2);
      if(diff < 0.1) is_one_peak = true;
      FloatType pi_180 = scitbx::constants::pi_180;
      cctbx_frac site_cart_h1 = site_cart_h1_fitted;
      cctbx_frac site_cart_h2 = site_cart_h2_fitted;
      for(FloatType x=0; x<360; x+=angular_shift) {
        for(FloatType y=0; y<360; y+=angular_shift) {
          for(FloatType z=0; z<360; z+=angular_shift) {
            FloatType x_ = x * pi_180;
            FloatType y_ = y * pi_180;
            FloatType z_ = z * pi_180;
            scitbx::af::tiny<double,3> c(std::cos(x_),std::cos(y_),std::cos(z_));
            scitbx::af::tiny<double,3> s(std::sin(x_),std::sin(y_),std::sin(z_));
            mat3<double> rot_mat = mat3<double>(
               c[0]*c[1]*c[2]-s[0]*s[2],
              -c[0]*c[1]*s[2]-s[0]*c[2],
               c[0]*s[1],
               s[0]*c[1]*c[2]+c[0]*s[2],
              -s[0]*c[1]*s[2]+c[0]*c[2],
               s[0]*s[1],
              -s[1]*c[2],
               s[1]*s[2],
               c[1]);
            cctbx_frac sites_cart_h1_new =
              (site_cart_h1 - origin_cart) * rot_mat + origin_cart;
            cctbx_frac sites_cart_h2_new =
              (site_cart_h2 - origin_cart) * rot_mat + origin_cart;
            cctbx_frac sites_frac_h1_new = unit_cell.fractionalize(
              sites_cart_h1_new);
            FloatType dist = std::sqrt(unit_cell.distance_sq(sites_frac_h1_new,
              site_frac_peak1));
            if(!is_one_peak) {
              cctbx_frac sites_frac_h2_new = unit_cell.fractionalize(
                sites_cart_h2_new);
              dist += std::sqrt(unit_cell.distance_sq(sites_frac_h2_new,site_frac_peak2));
            }
            if(dist < dist_best) {
              dist_best = dist;
              site_cart_o_fitted = origin_cart;
              site_cart_h1_fitted = sites_cart_h1_new;
              site_cart_h2_fitted = sites_cart_h2_new;
            }
          }
        }
      }
    }
};

}} // namespace mmtbx::utils

#endif // MMTBX_UTILS_H
