#ifndef CCTBX_XRAY_SCATTERING_DICTIONARY_H
#define CCTBX_XRAY_SCATTERING_DICTIONARY_H

#include <cctbx/xray/scatterer.h>
#include <cctbx/eltbx/caasf.h>
#include <map>

namespace cctbx { namespace xray {

  struct scatterer_group
  {
    scatterer_group() : coefficients(0) {}

    eltbx::caasf::custom coefficients;
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
      find_all_zero() const
      {
        af::shared<std::string> result;
        for(dict_type::const_iterator e=dict_.begin();e!=dict_.end();e++) {
          if (e->second.coefficients.all_zero()) {
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
        eltbx::caasf::custom const& coefficients)
      {
        CCTBX_ASSERT(dict_.find(label) != dict_.end());
        dict_[label].coefficients = coefficients;
      }

      void
      assign_from_table(std::string const& table)
      {
        CCTBX_ASSERT(table == "IT1992" || table == "WK1995");
        if (table == "IT1992") {
          for(dict_type::iterator e=dict_.begin();e!=dict_.end();e++) {
            if (!process_const_and_custom(e)) {
              eltbx::caasf::it1992 entry(e->first, 1);
              e->second.coefficients = eltbx::caasf::custom(
                entry.a(), entry.b(), entry.c());
            }
          }
        }
        else {
          for(dict_type::iterator e=dict_.begin();e!=dict_.end();e++) {
            if (!process_const_and_custom(e)) {
              eltbx::caasf::wk1995 entry(e->first, 1);
              e->second.coefficients = eltbx::caasf::custom(
                entry.a(), entry.b(), entry.c());
            }
          }
        }
      }

    protected:
      std::size_t n_scatterers_;
      dict_type dict_;

      bool
      process_const_and_custom(dict_type::iterator e)
      {
        if (e->first == "const") {
          e->second.coefficients = eltbx::caasf::custom(1);
          return true;
        }
        if (e->first == "custom") {
          return true;
        }
        return false;
      }
  };

}} // namespace cctbx::xray

#endif // CCTBX_XRAY_SCATTERING_DICTIONARY_H
