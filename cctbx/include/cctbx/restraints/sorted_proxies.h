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
        asu_mappings_(asu_mappings.get())
      {}

      //! Instance as passed to the constructor.
      boost::shared_ptr<direct_space_asu::asu_mappings<> > const&
      asu_mappings() const { return asu_mappings_owner_; }

      bool
      process(ProxyType const& proxy)
      {
        proxies.push_back(proxy);
        return false;
      }

      bool
      process(SymProxyType const& proxy)
      {
        if (asu_mappings_->is_direct_interaction(proxy.pair)) {
          if (proxy.pair.j_sym == 0 || proxy.pair.i_seq < proxy.pair.j_seq) {
            proxies.push_back(proxy.as_direct_proxy());
          }
          return false;
        }
        sym_proxies.push_back(proxy);
        return true;
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
  };

}} // namespace cctbx::restraints

#endif // CCTBX_RESTRAINTS_SORTED_PROXIES_H
