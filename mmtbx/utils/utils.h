#ifndef MMTBX_UTILS_H
#define MMTBX_UTILS_H

#include <scitbx/array_family/shared.h>
#include <mmtbx/error.h>
#include <cctbx/uctbx.h>

using namespace std;
namespace mmtbx { namespace utils {
namespace af=scitbx::af;
using scitbx::mat3;
using scitbx::vec3;
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

template <typename FloatType>
af::shared<std::size_t>
  select_water_by_distance(af::shared<vec3<FloatType> > const& sites_frac_all,
                           af::shared<std::string> const& element_symbols_all,
                           af::shared<std::size_t> const& water_selection_o,
                           FloatType const& dist_max,
                           FloatType const& dist_min_mac,
                           FloatType const& dist_min_sol,
                           cctbx::uctbx::unit_cell const& unit_cell)
{
  af::shared<std::size_t> result_selection;
  std::size_t closest_index=0;
  for(std::size_t i=0; i<water_selection_o.size(); i+=1) {
    std::size_t i_wat = water_selection_o[i];
    MMTBX_ASSERT(element_symbols_all[i_wat] != "H");
    MMTBX_ASSERT(element_symbols_all[i_wat] != "D");
    FloatType dist_closest = 999.;
    std::string closest_element;
    for(std::size_t j=0; j<sites_frac_all.size(); j+=1) {
      if(element_symbols_all[j]!="H"&&element_symbols_all[j]!="D"&&j!=i_wat) {
        FloatType dist = unit_cell.distance(
          cctbx::fractional<>(sites_frac_all[i_wat]),
          cctbx::fractional<>(sites_frac_all[j]));
        if(dist < dist_closest) {
          dist_closest = dist;
          closest_element = element_symbols_all[j];
          closest_index = j;
        }
      }
    }
    bool is_closest_index_water = false;
    for(std::size_t k=0; k<water_selection_o.size(); k+=1) {
      if(water_selection_o[k] == closest_index) {
        is_closest_index_water = true;
        break;
      }
    }
    double dist_min = dist_min_mac;
    if(is_closest_index_water) dist_min = dist_min_sol;
    if(dist_closest<=dist_max&&dist_closest>=dist_min&&(closest_element=="N"||
       closest_element=="O")) {
      result_selection.push_back(i_wat);
    }
  }
  return result_selection;
}

///
template <typename FloatType>
void correct_drifted_waters(af::ref<vec3<FloatType> > const& sites_frac_all,
                            af::const_ref<vec3<FloatType> > const& sites_frac_peaks,
                            af::const_ref<bool> const& water_selection,
                            cctbx::uctbx::unit_cell const& unit_cell)
{
  MMTBX_ASSERT(sites_frac_all.size() == water_selection.size());
  for(std::size_t i=0; i<sites_frac_all.size(); i+=1) {
    if(water_selection[i]) {
      FloatType dist_min = 999.;
      vec3<FloatType> closest_site;
      for(std::size_t j=0; j<sites_frac_peaks.size(); j+=1) {
        FloatType dist = unit_cell.distance(
          cctbx::fractional<>(sites_frac_all[i]),
          cctbx::fractional<>(sites_frac_peaks[j]));
        if(dist < dist_min) {
          dist_min = dist;
          closest_site = sites_frac_peaks[j];
        }
      }
      if(dist_min < 0.5 && dist_min > 0.1) {
        sites_frac_all[i] = closest_site;
      }
    }
  }
}
///


}} // namespace mmtbx::utils

#endif // MMTBX_UTILS_H
