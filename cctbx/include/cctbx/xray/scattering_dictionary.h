#ifndef CCTBX_XRAY_SCATTERING_DICTIONARY_H
#define CCTBX_XRAY_SCATTERING_DICTIONARY_H

#include <cctbx/xray/scatterer.h>
#include <cctbx/eltbx/xray_scattering.h>
#include <map>

namespace cctbx { namespace xray {

  struct scatterer_group
  {
    scatterer_group() : gaussian(0) {}

    eltbx::xray_scattering::gaussian gaussian;
    af::shared<std::size_t> member_indices;
  };

  class scattering_dictionary
  {
    public:
      typedef std::map<std::string, scatterer_group> dict_type;

      scattering_dictionary()
      :
        n_scatterers_(0)
      {}

      template <typename XrayScattererType>
      scattering_dictionary(
        af::const_ref<XrayScattererType> const& xray_scatterers)
      :
        n_scatterers_(xray_scatterers.size())
      {
        for(std::size_t i=0;i<xray_scatterers.size();i++) {
          dict_[xray_scatterers[i].scattering_type]
            .member_indices.push_back(i);
        }
      }

      std::size_t
      n_scatterers() const { return n_scatterers_; }

      dict_type const&
      dict() const { return dict_; }

      dict_type&
      dict() { return dict_; }

      af::shared<std::string>
      find_undefined() const
      {
        af::shared<std::string> result;
        for(dict_type::const_iterator e=dict_.begin();e!=dict_.end();e++) {
          if (e->second.gaussian.all_zero()) {
            result.push_back(e->first);
          }
        }
        return result;
      }

      scatterer_group const&
      lookup(std::string const& label) const
      {
        dict_type::const_iterator e = dict_.find(label);
        if (e == dict_.end()) {
          throw error("Label not in scattering dictionary: " + label);
        }
        return e->second;
      }

      af::shared<std::size_t>
      scatterer_permutation() const
      {
        af::shared<std::size_t> result((af::reserve(n_scatterers_)));
        for(dict_type::const_iterator e=dict_.begin();e!=dict_.end();e++) {
          af::const_ref<std::size_t>
            member_indices = e->second.member_indices.const_ref();
          for(std::size_t mi=0;mi<member_indices.size();mi++) {
            result.push_back(member_indices[mi]);
          }
        }
        return result;
      }

      void
      assign(
        std::string const& label,
        eltbx::xray_scattering::gaussian const& gaussian)
      {
        CCTBX_ASSERT(dict_.find(label) != dict_.end());
        dict_[label].gaussian = gaussian;
      }

      void
      assign_from_table(std::string const& table)
      {
        CCTBX_ASSERT(table == "IT1992" || table == "WK1995");
        if (table == "IT1992") {
          for(dict_type::iterator e=dict_.begin();e!=dict_.end();e++) {
            if (!e->second.gaussian.all_zero()) continue;
            e->second.gaussian
              = eltbx::xray_scattering::it1992(e->first, 1).fetch();
          }
        }
        else {
          for(dict_type::iterator e=dict_.begin();e!=dict_.end();e++) {
            if (!e->second.gaussian.all_zero()) continue;
            e->second.gaussian
              = eltbx::xray_scattering::wk1995(e->first, 1).fetch();
          }
        }
      }

    protected:
      std::size_t n_scatterers_;
      dict_type dict_;
  };

}} // namespace cctbx::xray

#endif // CCTBX_XRAY_SCATTERING_DICTIONARY_H
