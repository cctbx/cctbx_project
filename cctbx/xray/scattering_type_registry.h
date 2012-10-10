#ifndef CCTBX_XRAY_SCATTERING_TYPE_REGISTRY_H
#define CCTBX_XRAY_SCATTERING_TYPE_REGISTRY_H

#include <cctbx/xray/scatterer.h>
#include <cctbx/eltbx/xray_scattering.h>
#include <cctbx/eltbx/electron_scattering.h>
#include <cctbx/eltbx/neutron.h>
#include <boost/optional.hpp>
#include <map>

namespace cctbx { namespace xray {

  class scattering_type_registry
  {
    public:
      typedef std::map<std::string, std::size_t> type_index_pairs_t;
      typedef eltbx::xray_scattering::gaussian gaussian_t;
      typedef af::shared<boost::optional<eltbx::xray_scattering::gaussian> >
              unique_gaussians_t;
      typedef af::shared<std::size_t> unique_counts_t;
      type_index_pairs_t type_index_pairs;
      unique_gaussians_t unique_gaussians;
      unique_counts_t unique_counts;

      static
      std::runtime_error not_in_registry(std::string const &scattering_type) {
        return std::runtime_error("scattering_type \""
                                  + scattering_type
                                  + "\" not in scattering_type_registry.");
      }

      static
      std::runtime_error
      form_factor_not_defined(std::string const &scattering_type) {
        return std::runtime_error(
          "gaussian not defined for scattering_type \""
          + scattering_type
          + "\".");
      }

      scattering_type_registry() {}

      std::size_t
      size() const
      {
        CCTBX_ASSERT(unique_gaussians.size() == type_index_pairs.size());
        CCTBX_ASSERT(unique_counts.size() == type_index_pairs.size());
        return type_index_pairs.size();
      }

      bool
      has_key(std::string const& scattering_type) const
      {
        return (   type_index_pairs.find(scattering_type)
                != type_index_pairs.end());
      }

      std::size_t
      process(std::string const& scattering_type)
      {
        type_index_pairs_t::const_iterator
          pair = type_index_pairs.find(scattering_type);
        if (pair != type_index_pairs.end()) {
          unique_counts[pair->second]++;
          return pair->second;
        }
        std::size_t index = unique_gaussians.size();
        type_index_pairs[scattering_type] = index;
        unique_gaussians.push_back(boost::optional<gaussian_t>());
        unique_counts.push_back(1);
        return index;
      }

      template <typename XrayScattererType>
      af::shared<std::size_t>
      process(af::const_ref<XrayScattererType> const& scatterers)
      {
        af::shared<std::size_t> result(
          scatterers.size(), af::init_functor_null<std::size_t>());
        for(std::size_t i=0;i<scatterers.size();i++) {
          result[i] = process(scatterers[i].scattering_type);
        }
        return result;
      }

      std::size_t
      unique_index(std::string const& scattering_type) const
      {
        type_index_pairs_t::const_iterator
          pair = type_index_pairs.find(scattering_type);
        if (pair != type_index_pairs.end()) return pair->second;
        throw not_in_registry(scattering_type);
      }

      template <typename XrayScattererType>
      af::shared<std::size_t>
      unique_indices(af::const_ref<XrayScattererType> const& scatterers) const
      {
        af::shared<std::size_t> result(
          scatterers.size(), af::init_functor_null<std::size_t>());
        for(std::size_t i=0;i<scatterers.size();i++) {
          result[i] = unique_index(scatterers[i].scattering_type);
        }
        return result;
      }

      template <typename XrayScattererType>
      af::shared<double>
      occupancy_sums(af::const_ref<XrayScattererType> const& scatterers) const
      {
        af::shared<double> result(unique_counts.size(), 0);
        for(std::size_t i=0;i<scatterers.size();i++) {
          XrayScattererType const& sc = scatterers[i];
          result[unique_index(sc.scattering_type)] += sc.occupancy;
        }
        return result;
      }

      template <typename XrayScattererType>
      af::shared<double>
      unit_cell_occupancy_sums(af::const_ref<XrayScattererType> const&
                               scatterers) const
      {
        af::shared<double> result(unique_counts.size(), 0);
        for(std::size_t i=0;i<scatterers.size();i++) {
          XrayScattererType const& sc = scatterers[i];
          result[unique_index(sc.scattering_type)
                 ] += sc.occupancy*sc.multiplicity();
        }
        return result;
      }

      boost::optional<eltbx::xray_scattering::gaussian> const&
      gaussian(std::string const& scattering_type) const
      {
        return unique_gaussians[unique_index(scattering_type)];
      }

      eltbx::xray_scattering::gaussian const&
      gaussian_not_optional(std::string const& scattering_type) const
      {
        boost::optional<gaussian_t> const& result = gaussian(scattering_type);
        if (!result) {
          throw form_factor_not_defined(scattering_type);
        }
        return *result;
      }

      af::shared<std::string>
      unassigned_types() const
      {
        af::shared<std::string> result;
        af::const_ref<boost::optional<gaussian_t> >
          ugs = unique_gaussians.const_ref();
        for(type_index_pairs_t::const_iterator
              pair=type_index_pairs.begin();
              pair!=type_index_pairs.end();
              pair++) {
          std::size_t ui = pair->second;
          if (!ugs[ui]) result.push_back(pair->first);
        }
        return result;
      }

      bool
      assign(
        std::string const& scattering_type,
        boost::optional<scitbx::math::gaussian::sum<double> > const& gaussian)
      {
        std::size_t ui = unique_index(scattering_type);
        bool result = !unique_gaussians[ui];
        if (!gaussian) unique_gaussians[ui] = boost::optional<gaussian_t>();
        else           unique_gaussians[ui] = gaussian_t(*gaussian);
        return result;
      }

      void
      assign_from_table(std::string const& table, bool exact = true)
      {
        CCTBX_ASSERT(table == "IT1992" || table == "WK1995" ||
                                 table == "PENG1996" || table == "NEUTRON1992");
        bool has_assigned = false;
        af::ref<boost::optional<gaussian_t> > ugs = unique_gaussians.ref();
        if (table == "IT1992") {
          for(type_index_pairs_t::const_iterator
                pair=type_index_pairs.begin();
                pair!=type_index_pairs.end();
                pair++) {
            std::size_t ui = pair->second;
            if (ugs[ui]) continue;
            ugs[ui] = eltbx::xray_scattering::it1992(
              pair->first, exact).fetch();
            has_assigned = true;
          }
        }
        else if (table == "PENG1996") {
          for(type_index_pairs_t::const_iterator
                pair=type_index_pairs.begin();
                pair!=type_index_pairs.end();
                pair++) {
            std::size_t ui = pair->second;
            if (ugs[ui]) continue;
            ugs[ui] = eltbx::electron_scattering::peng1996(
              pair->first, exact).fetch();
            has_assigned = true;
          }
        }
        else if (table == "NEUTRON1992") {
          for(type_index_pairs_t::const_iterator
                pair=type_index_pairs.begin();
                pair!=type_index_pairs.end();
                pair++) {
            std::size_t ui = pair->second;
            if (ugs[ui]) continue;
            ugs[ui] = eltbx::neutron::neutron_news_1992_table(
              pair->first, exact).fetch();
            has_assigned = true;
          }
        }
        else {
          for(type_index_pairs_t::const_iterator
                pair=type_index_pairs.begin();
                pair!=type_index_pairs.end();
                pair++) {
            std::size_t ui = pair->second;
            if (ugs[ui]) continue;
            ugs[ui] = eltbx::xray_scattering::wk1995(
              pair->first, exact).fetch();
            has_assigned = true;
          }
        }
        if( has_assigned )
          this->set_last_table(table);
      }

      std::string
      type_given_unique_index(std::size_t unique_index) const
      {
        for(type_index_pairs_t::const_iterator
              pair=type_index_pairs.begin();
              pair!=type_index_pairs.end();
              pair++) {
          if (pair->second == unique_index) return pair->first;
        }
        throw std::runtime_error("unique_index out of range.");
      }

      af::shared<double>
      unique_form_factors_at_d_star_sq(double d_star_sq) const
      {
        af::const_ref<boost::optional<gaussian_t> >
          ugs = unique_gaussians.const_ref();
        af::shared<double> result(ugs.size(), af::init_functor_null<double>());
        double x_sq = d_star_sq / 4;
        for(std::size_t i=0;i<ugs.size();i++) {
          if (!ugs[i]) {
            throw form_factor_not_defined(type_given_unique_index(i));
          }
          result[i] = ugs[i]->at_x_sq(x_sq);
        }
        return result;
      }

      af::shared<double>
      dilated_form_factors_at_d_star_sq(
        double d_star_sq,
        af::const_ref<double> const &dilation_coefficients,
        unique_counts_t unique_indices)
      {
        CCTBX_ASSERT(dilation_coefficients.size() == unique_indices.size());
        af::shared<double> result(dilation_coefficients.size());
        af::const_ref< boost::optional<gaussian_t> >
          ugs = unique_gaussians.const_ref();
        for(std::size_t i=0; i < dilation_coefficients.size(); ++i) {
          std::size_t j = unique_indices[i];
          boost::optional<gaussian_t> const &ff = ugs[j];
          if (!ff) throw form_factor_not_defined(type_given_unique_index(j));
          result[i] = ff->at_d_star_sq(d_star_sq/dilation_coefficients[i]);
        }
        return result;
      }

      void set_last_table(const std::string &table)
      {
        last_table_ = table;
      }

      const std::string& last_table() const {return last_table_;}

    private:
      std::string last_table_;
  };

}} // namespace cctbx::xray

#endif // CCTBX_XRAY_SCATTERING_TYPE_REGISTRY_H
