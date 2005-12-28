#ifndef CCTBX_GEOMETRY_RESTRAINTS_ASU_CACHE_H
#define CCTBX_GEOMETRY_RESTRAINTS_ASU_CACHE_H

#include <cctbx/crystal/direct_space_asu.h>

namespace cctbx { namespace geometry_restraints {

  namespace direct_space_asu = crystal::direct_space_asu;

  //! Asymmetric unit cache to facilitate speed optimizations.
  /*! Not available in Python.
   */
  template <typename FloatType=double, typename IntShiftType=int>
  class asu_cache
  {
    public:
      //! Convenience typedef.
      typedef direct_space_asu::asu_mappings<FloatType, IntShiftType>
        asu_mappings_t;

      //! Default constructor. Some data members are not initialized!
      asu_cache() {}

      //! Support for bond_residual_sum and nonbonded_residual_sum functions.
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
        scitbx::vec3<FloatType>* sites_ptr = (
          sites_memory_.size() ? &*sites_memory_.begin() : 0);
        std::size_t sum_n_sym = 0;
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
            sum_n_sym += n_sym;
          }
        }
        CCTBX_ASSERT(sum_n_sym <= sites_memory_.size());
        if (allocate_gradients) {
          gradients.resize(n_sites, scitbx::vec3<FloatType>(0,0,0));
        }
      }

      //! Support for bond_residual_sum and nonbonded_residual_sum functions.
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

      /*! \brief Array of cartesian coordinates inside the asymmetric unit
          and the buffer region.
       */
      std::vector<scitbx::vec3<FloatType>*> sites;
      //! Array of gradients.
      std::vector<scitbx::vec3<FloatType> > gradients;

    protected:
      std::vector<scitbx::vec3<FloatType> > sites_memory_;
      af::const_ref<typename asu_mappings_t::array_of_mappings_for_one_site>
        mappings_;
  };

}} // namespace cctbx::geometry_restraints

#endif // CCTBX_GEOMETRY_RESTRAINTS_ASU_CACHE_H
