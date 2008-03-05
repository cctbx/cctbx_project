#ifndef MMTBX_UTILS_H
#define MMTBX_UTILS_H

#include <scitbx/array_family/shared.h>
#include <mmtbx/error.h>
#include <cctbx/uctbx.h>

using namespace std;
namespace mmtbx { namespace utils {
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
    FloatType dist_best_sq;

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
    dist_best_sq(1.e+9)
    {
      CCTBX_ASSERT(angular_shift > 0 && angular_shift < 360);
      bool is_one_peak = false;
      FloatType diff = unit_cell.distance_sq(site_frac_peak1, site_frac_peak2);
      if(diff < 0.1) is_one_peak = true;
      FloatType pi_180 = scitbx::constants::pi_180;
      cctbx_frac site_cart_h1 = site_cart_h1_fitted;
      cctbx_frac site_cart_h2 = site_cart_h2_fitted;
      for(FloatType x=0; x<360; x+=angular_shift) {
        FloatType x_ = x * pi_180;
        FloatType cos_x = std::cos(x_);
        FloatType sin_x = std::sin(x_);
        for(FloatType y=0; y<360; y+=angular_shift) {
          FloatType y_ = y * pi_180;
          FloatType cos_y = std::cos(y_);
          FloatType sin_y = std::sin(y_);
          for(FloatType z=0; z<360; z+=angular_shift) {
            FloatType z_ = z * pi_180;
            FloatType cos_z = std::cos(z_);
            FloatType sin_z = std::sin(z_);
            mat3<double> rot_mat = mat3<double>(
               cos_x*cos_y*cos_z-sin_x*sin_z,
              -cos_x*cos_y*sin_z-sin_x*cos_z,
               cos_x*sin_y,
               sin_x*cos_y*cos_z+cos_x*sin_z,
              -sin_x*cos_y*sin_z+cos_x*cos_z,
               sin_x*sin_y,
              -sin_y*cos_z,
               sin_y*sin_z,
               cos_y);
            cctbx_frac sites_cart_h1_new =
              (site_cart_h1 - origin_cart) * rot_mat + origin_cart;
            cctbx_frac sites_cart_h2_new =
              (site_cart_h2 - origin_cart) * rot_mat + origin_cart;
            cctbx_frac sites_frac_h1_new = unit_cell.fractionalize(
              sites_cart_h1_new);
            FloatType dist = unit_cell.distance_sq(sites_frac_h1_new,
              site_frac_peak1);
            if(!is_one_peak) {
              cctbx_frac sites_frac_h2_new = unit_cell.fractionalize(
                sites_cart_h2_new);
              dist += unit_cell.distance_sq(sites_frac_h2_new,site_frac_peak2);
            }
            if(dist < dist_best_sq) {
              dist_best_sq = dist;
              site_cart_o_fitted = origin_cart;
              site_cart_h1_fitted = sites_cart_h1_new;
              site_cart_h2_fitted = sites_cart_h2_new;
            }
          }
        }
      }
    }

    double dist_best() { return std::sqrt(dist_best_sq); }
};

}} // namespace mmtbx::utils

#endif // MMTBX_UTILS_H
