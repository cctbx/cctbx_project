#ifndef CCTBX_RESTRAINTS_SORTED_PROXIES_H
#define CCTBX_RESTRAINTS_SORTED_PROXIES_H

#include <cctbx/crystal/direct_space_asu.h>

namespace cctbx { namespace restraints {

  template <typename ProxyType,
            typename SymProxyType>
  class sorted_proxies
  {
    public:
      sorted_proxies() {}

      sorted_proxies(
        boost::shared_ptr<
          direct_space_asu::asu_mappings<> > const& asu_mappings)
      :
        asu_mappings_owner_(asu_mappings),
        asu_mappings_(asu_mappings.get()),
        sym_active_flags(asu_mappings->mappings_const_ref().size(), false)
      {}

      //! Instance as passed to the constructor.
      boost::shared_ptr<direct_space_asu::asu_mappings<> > const&
      asu_mappings() const { return asu_mappings_owner_; }

      bool
      process(ProxyType const& proxy)
      {
        CCTBX_ASSERT(proxy.i_seqs[0] < proxy.i_seqs[1]);
        proxies.push_back(proxy);
        return false;
      }

      bool
      process(SymProxyType const& proxy)
      {
        int type_id = asu_mappings_->interaction_type_id(proxy.pair);
        if (type_id < 0) return false;
        if (type_id > 0) {
          proxies.push_back(proxy.as_simple_proxy());
          return false;
        }
        push_back(proxy);
        return true;
      }

      void
      push_back(SymProxyType const& proxy)
      {
        sym_proxies.push_back(proxy);
        sym_active_flags[proxy.pair.i_seq] = true;
        sym_active_flags[proxy.pair.j_seq] = true;
      }

      std::size_t
      n_total() const
      {
        return proxies.size() + sym_proxies.size();
      }

    protected:
      boost::shared_ptr<direct_space_asu::asu_mappings<> > asu_mappings_owner_;
      const direct_space_asu::asu_mappings<>* asu_mappings_;

    public:
      af::shared<ProxyType> proxies;
      af::shared<SymProxyType> sym_proxies;
      std::vector<bool> sym_active_flags;
  };

}} // namespace cctbx::restraints

#endif // CCTBX_RESTRAINTS_SORTED_PROXIES_H
