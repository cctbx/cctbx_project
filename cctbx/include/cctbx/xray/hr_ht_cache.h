#ifndef CCTBX_XRAY_HR_HT_CACHE_H
#define CCTBX_XRAY_HR_HT_CACHE_H

#include <cctbx/sgtbx/miller_ops.h>

namespace cctbx { namespace xray { namespace structure_factors {

  template <typename FloatType>
  struct hr_ht_group
  {
    hr_ht_group(
      miller::index<> const& hr_,
      FloatType const& ht_)
    :
      hr(hr_),
      ht(ht_)
    {}

    miller::index<> hr;
    FloatType ht;
  };

  template <typename FloatType>
  struct hr_ht_cache
  {
    typedef FloatType f_t;

    hr_ht_cache(
      sgtbx::space_group const& space_group,
      miller::index<> const& h)
    :
      is_centric(space_group.is_centric())
    {
      f_t t_den = static_cast<f_t>(space_group.t_den());
      if (!is_centric) {
        h_inv_t = -1;
        is_origin_centric = false;
      }
      else {
        h_inv_t = static_cast<f_t>(h * space_group.inv_t()) / t_den;
        is_origin_centric = (h_inv_t == 0);
      }
      for(std::size_t i_smx=0;i_smx<space_group.n_smx();i_smx++) {
        sgtbx::rt_mx const& s = space_group.smx(i_smx);
        groups.push_back(hr_ht_group<f_t>(
          h * s.r(),
          static_cast<f_t>(h * s.t()) / t_den));
      }
    }

    bool is_centric;
    bool is_origin_centric;
    f_t h_inv_t;
    af::small<hr_ht_group<f_t>, sgtbx::n_max_repr_rot_mx> groups;
  };

}}} // namespace cctbx::xray::structure_factors

#endif // CCTBX_XRAY_HR_HT_CACHE_H
