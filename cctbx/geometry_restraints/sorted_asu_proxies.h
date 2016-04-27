#ifndef CCTBX_GEOMETRY_RESTRAINTS_SORTED_ASU_PROXIES_H
#define CCTBX_GEOMETRY_RESTRAINTS_SORTED_ASU_PROXIES_H

#include <cctbx/crystal/direct_space_asu.h>

namespace cctbx { namespace geometry_restraints {

  //! Convenience typedef.
  typedef direct_space_asu::asu_mapping_index      asu_mapping_index;
  //! Convenience typedef.
  typedef direct_space_asu::asu_mapping_index_pair asu_mapping_index_pair;
  //! Convenience typedef.
  typedef direct_space_asu::asu_mappings<>         asu_mappings;

  //! Managed group of simple proxies and asu proxies.
  /*! Sorting of the proxies enables optimizations for speed.

      See also:
        bond_sorted_asu_proxies,
        nonbonded_sorted_asu_proxies
   */
  template <typename SimpleProxyType,
            typename AsuProxyType>
  class sorted_asu_proxies
  {
    public:
      //! Default constructor. Some data members are not initialized!
      sorted_asu_proxies()
      :
        asu_mappings_(0)
      {}

      //! Initialization with asu_mappings.
      sorted_asu_proxies(
        boost::shared_ptr<
          direct_space_asu::asu_mappings<> > const& asu_mappings)
      :
        asu_mappings_owner_(asu_mappings),
        asu_mappings_(asu_mappings.get())
      {
        if (asu_mappings_ != 0) {
          asu_active_flags.resize(
            asu_mappings_->mappings_const_ref().size(), false);
        }
      }

      //! Instance as passed to the constructor.
      boost::shared_ptr<direct_space_asu::asu_mappings<> > const&
      asu_mappings() const
      {
        CCTBX_ASSERT(asu_mappings_ != 0);
        return asu_mappings_owner_;
      }

      //! Appends proxy to simple array.
      void
      process(SimpleProxyType const& proxy)
      {
        simple.push_back(proxy);
      }

      //! Appends all proxies to simple array.
      void
      process(af::const_ref<SimpleProxyType> const& proxies)
      {
        for(std::size_t i=0;i<proxies.size();i++) process(proxies[i]);
      }

      /*! \brief Appends proxy to the simple array if possible, the
          asu array otherwise.
       */
      /*! See also:
            asu_mappings::is_simple_interaction,
            bond_asu_proxy::as_simple_proxy,
            nonbonded_asu_proxy::as_simple_proxy
       */
      void
      process(AsuProxyType const& proxy, bool sym_excl_flag=false)
      {
        CCTBX_ASSERT(asu_mappings_ != 0 && proxy.is_active());
        if (asu_mappings_->is_simple_interaction(proxy)) {
          if (proxy.i_seq < proxy.j_seq) {
            simple.push_back(proxy.as_simple_proxy());
          }
        }
        else if (!sym_excl_flag) {
          push_back(proxy);
        }
      }

      //! Calls process() for each proxy in the proxies array.
      void
      process(af::const_ref<AsuProxyType> const& proxies)
      {
        for(std::size_t i=0;i<proxies.size();i++) process(proxies[i]);
      }

      //! Unconditionally appends proxy to the asu array (for debugging).
      void
      push_back(AsuProxyType const& proxy)
      {
        CCTBX_ASSERT(asu_mappings_ != 0);
        CCTBX_ASSERT(proxy.i_seq < asu_active_flags.size());
        CCTBX_ASSERT(proxy.j_seq < asu_active_flags.size());
        asu.push_back(proxy);
        asu_active_flags[proxy.i_seq] = true;
        asu_active_flags[proxy.j_seq] = true;
      }

      //! Unconditionally appends all proxies to the asu array (for debugging).
      void
      push_back(af::const_ref<AsuProxyType> const& proxies)
      {
        for(std::size_t i=0;i<proxies.size();i++) push_back(proxies[i]);
      }

      //! Shorthand for: simple.size() + asu.size()
      std::size_t
      n_total() const { return simple.size() + asu.size(); }

//    protected:
    public:
      boost::shared_ptr<direct_space_asu::asu_mappings<> > asu_mappings_owner_;
      const direct_space_asu::asu_mappings<>* asu_mappings_;

    public:
      //! Array of simple proxies.
      af::shared<SimpleProxyType> simple;
      //! Array of asu proxies.
      af::shared<AsuProxyType> asu;
      /*! \brief Array of size asu_mappings().mappings().size() indicating
          if a site is involved in one or more asu proxies.
       */
      std::vector<bool> asu_active_flags;
  };

}} // namespace cctbx::geometry_restraints

#endif // CCTBX_GEOMETRY_RESTRAINTS_SORTED_ASU_PROXIES_H
