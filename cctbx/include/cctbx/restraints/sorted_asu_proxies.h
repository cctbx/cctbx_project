#ifndef CCTBX_RESTRAINTS_SORTED_ASU_PROXIES_H
#define CCTBX_RESTRAINTS_SORTED_ASU_PROXIES_H

#include <cctbx/crystal/direct_space_asu.h>

namespace cctbx { namespace restraints {

  typedef direct_space_asu::asu_mapping_index      asu_mapping_index;
  typedef direct_space_asu::asu_mapping_index_pair asu_mapping_index_pair;
  typedef direct_space_asu::asu_mappings<>         asu_mappings;

  template <typename SimpleProxyType,
            typename SymProxyType>
  class sorted_asu_proxies
  {
    public:
      sorted_asu_proxies() {}

      sorted_asu_proxies(
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
      process(SimpleProxyType const& proxy)
      {
        simple.push_back(proxy);
        return false;
      }

      bool
      process(SymProxyType const& proxy)
      {
        if (asu_mappings_->is_simple_interaction(proxy)) {
          simple.push_back(proxy.as_simple_proxy());
          return false;
        }
        push_back(proxy);
        return true;
      }

      void
      push_back(SymProxyType const& proxy)
      {
        sym.push_back(proxy);
        sym_active_flags[proxy.i_seq] = true;
        sym_active_flags[proxy.j_seq] = true;
      }

      std::size_t
      n_total() const { return simple.size() + sym.size(); }

    protected:
      boost::shared_ptr<direct_space_asu::asu_mappings<> > asu_mappings_owner_;
      const direct_space_asu::asu_mappings<>* asu_mappings_;

    public:
      af::shared<SimpleProxyType> simple;
      af::shared<SymProxyType> sym;
      std::vector<bool> sym_active_flags;
  };

}} // namespace cctbx::restraints

#endif // CCTBX_RESTRAINTS_SORTED_ASU_PROXIES_H
