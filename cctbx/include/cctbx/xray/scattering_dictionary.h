#ifndef CCTBX_XRAY_SCATTERING_DICTIONARY_H
#define CCTBX_XRAY_SCATTERING_DICTIONARY_H

#include <cctbx/xray/scatterer.h>
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

      scattering_dictionary() {}

      template <typename XrayScattererType>
      scattering_dictionary(
        af::const_ref<XrayScattererType> const& xray_scatterers)
      {
        for(std::size_t i=0;i<xray_scatterers.size();i++) {
          dict_[xray_scatterers[i].caasf.label()].member_indices.push_back(i);
        }
      }

      std::size_t
      size() const { return dict_.size(); }

      dict_type const&
      dict() const { return dict_; }

      dict_type&
      dict() { return dict_; }

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
      dict_type dict_;

      bool
      process_const_and_custom(dict_type::iterator e)
      {
        if (e->first == "const") {
          e->second.coefficients = eltbx::caasf::custom(1);
          return true;
        }
        if (std::strncmp(e->first.c_str(), "custom", 6) == 0) {
          return true;
        }
        return false;
      }
  };

}} // namespace cctbx::xray

#endif // CCTBX_XRAY_SCATTERING_DICTIONARY_H
