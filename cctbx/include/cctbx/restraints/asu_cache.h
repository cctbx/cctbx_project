#ifndef CCTBX_RESTRAINTS_ASU_CACHE_H
#define CCTBX_RESTRAINTS_ASU_CACHE_H

#include <cctbx/crystal/direct_space_asu.h>

namespace cctbx { namespace restraints {

  namespace direct_space_asu = crystal::direct_space_asu;

  template <typename FloatType=double, typename IntShiftType=int>
  class asu_cache
  {
    public:
      typedef direct_space_asu::asu_mappings<FloatType, IntShiftType>
        asu_mappings_t;

      asu_cache() {}

      asu_cache(
        af::const_ref<scitbx::vec3<FloatType> > const& moved_sites_cart,
        asu_mappings_t const& asu_mappings,
        std::vector<bool> const& sym_active_flags,
        bool allocate_gradients)
      {
        std::size_t n_sites = moved_sites_cart.size();
        mappings_ = asu_mappings.mappings_const_ref();
        CCTBX_ASSERT(mappings_.size() == n_sites);
        sites_memory_.resize(asu_mappings.n_sites_in_asu_and_buffer());
        sites.resize(n_sites, 0);
        scitbx::vec3<FloatType>* sites_ptr = &*sites_memory_.begin();
        for(std::size_t i_seq=0;i_seq<n_sites;i_seq++) {
          if (!sym_active_flags[i_seq]) {
            sites[i_seq] = 0;
          }
          else {
            sites[i_seq] = sites_ptr;
            std::size_t n_sym = mappings_[i_seq].size();
            for(std::size_t i_sym=0;i_sym<n_sym;i_sym++) {
              *sites_ptr++ = asu_mappings.map_moved_site_to_asu(
                moved_sites_cart[i_seq], i_seq, i_sym);
            }
          }
        }
        CCTBX_ASSERT(sites_ptr <= &*sites_memory_.end());
        if (allocate_gradients) {
          gradients.resize(n_sites, scitbx::vec3<FloatType>(0,0,0));
        }
      }

      void
      add_gradients(
        af::ref<scitbx::vec3<double> > const& gradient_array,
        direct_space_asu::asu_mappings<> const& asu_mappings) const
      {
        for(std::size_t i_seq=0;i_seq<gradient_array.size();i_seq++) {
          gradient_array[i_seq] += asu_mappings.r_inv_cart(i_seq, 0)
                                 * gradients[i_seq];
        }
      }

      std::vector<scitbx::vec3<FloatType>*> sites;
      std::vector<scitbx::vec3<FloatType> > gradients;

    protected:
      std::vector<scitbx::vec3<FloatType> > sites_memory_;
      af::const_ref<typename asu_mappings_t::array_of_mappings_for_one_site>
        mappings_;
  };

}} // namespace cctbx::restraints

#endif // CCTBX_RESTRAINTS_ASU_CACHE_H
