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
      sorted_asu_proxies()
      :
        asu_mappings_(0)
      {}

      sorted_asu_proxies(
        boost::shared_ptr<
          direct_space_asu::asu_mappings<> > const& asu_mappings)
      :
        asu_mappings_owner_(asu_mappings),
        asu_mappings_(asu_mappings.get())
      {
        if (asu_mappings_ != 0) {
          sym_active_flags.resize(
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

      bool
      process(SimpleProxyType const& proxy)
      {
        simple.push_back(proxy);
        return false;
      }

      void
      process(af::const_ref<SimpleProxyType> const& proxies)
      {
        for(std::size_t i=0;i<proxies.size();i++) process(proxies[i]);
      }

      bool
      process(SymProxyType const& proxy)
      {
        CCTBX_ASSERT(asu_mappings_ != 0 && proxy.is_active());
        if (asu_mappings_->is_simple_interaction(proxy)) {
          if (proxy.i_seq < proxy.j_seq) {
            simple.push_back(proxy.as_simple_proxy());
          }
          return false;
        }
        push_back(proxy);
        return true;
      }

      void
      process(af::const_ref<SymProxyType> const& proxies)
      {
        for(std::size_t i=0;i<proxies.size();i++) process(proxies[i]);
      }

      void
      push_back(SymProxyType const& proxy)
      {
        CCTBX_ASSERT(asu_mappings_ != 0);
        CCTBX_ASSERT(proxy.i_seq < sym_active_flags.size());
        CCTBX_ASSERT(proxy.j_seq < sym_active_flags.size());
        sym.push_back(proxy);
        sym_active_flags[proxy.i_seq] = true;
        sym_active_flags[proxy.j_seq] = true;
      }

      void
      push_back(af::const_ref<SymProxyType> const& proxies)
      {
        for(std::size_t i=0;i<proxies.size();i++) push_back(proxies[i]);
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
