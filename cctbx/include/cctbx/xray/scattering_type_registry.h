#ifndef CCTBX_XRAY_SCATTERING_TYPE_REGISTRY_H
#define CCTBX_XRAY_SCATTERING_TYPE_REGISTRY_H

#include <cctbx/xray/scatterer.h>
#include <cctbx/eltbx/xray_scattering.h>
#include <boost/optional.hpp>
#include <map>

namespace cctbx { namespace xray {

  template<class FormFactorType>
  struct form_factor_table_policy;

  template<class FormFactorType>
  class generic_scattering_type_registry
  {
    public:
      typedef std::map<std::string, std::size_t> type_index_pairs_t;
      typedef FormFactorType form_factor_t;
      typedef af::shared< boost::optional<form_factor_t> >
              unique_form_factors_t;
      typedef af::shared<std::size_t> unique_counts_t;
      type_index_pairs_t type_index_pairs;
      unique_form_factors_t unique_form_factors;
      unique_counts_t unique_counts;

      generic_scattering_type_registry() {}

      std::size_t
      size() const
      {
        CCTBX_ASSERT(unique_form_factors.size() == type_index_pairs.size());
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
        std::size_t index = unique_form_factors.size();
        type_index_pairs[scattering_type] = index;
        unique_form_factors.push_back( boost::optional<form_factor_t>());
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
        throw std::runtime_error(
          "scattering_type \""
          + scattering_type
          + "\" not in scattering_type_registry.");
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

      boost::optional<form_factor_t> const&
      form_factor(std::string const& scattering_type) const
      {
        return unique_form_factors[unique_index(scattering_type)];
      }

      form_factor_t const&
      form_factor_not_optional(std::string const& scattering_type) const
      {
         boost::optional<form_factor_t> const&
           result = form_factor(scattering_type);
        if (!result) {
          throw std::runtime_error(
            "gaussian not defined for scattering_type \""
            + scattering_type
            + "\".");
        }
        return *result;
      }

      boost::optional<eltbx::xray_scattering::gaussian> const &
      gaussian(std::string const &scattering_type) const
      {
        return form_factor(scattering_type);
      }

      eltbx::xray_scattering::gaussian const &
      gaussian_not_optional(std::string const &scattering_type) const
      {
        return form_factor_not_optional(scattering_type);
      }

      af::shared<std::string>
      unassigned_types() const
      {
        af::shared<std::string> result;
        af::const_ref< boost::optional<form_factor_t> >
          uffs = unique_form_factors.const_ref();
        for(type_index_pairs_t::const_iterator
              pair=type_index_pairs.begin();
              pair!=type_index_pairs.end();
              pair++) {
          std::size_t ui = pair->second;
          if (!uffs[ui]) result.push_back(pair->first);
        }
        return result;
      }

      template<class CompatibleFormFactorType>
      bool
      assign(
        std::string const& scattering_type,
        boost::optional<CompatibleFormFactorType> const& form_factor)
      {
        std::size_t ui = unique_index(scattering_type);
        bool result = !unique_form_factors[ui];
        unique_form_factors[ui] = form_factor ? form_factor_t(*form_factor)
                                              : boost::optional<form_factor_t>();
        return result;
      }

      void
      assign_from_table(std::string const& table)
      {
        form_factor_table_policy<form_factor_t>::assign_from(
          table,
          unique_form_factors.ref(),
          type_index_pairs.begin(),
          type_index_pairs.end());
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
        af::const_ref< boost::optional<form_factor_t> >
          uffs = unique_form_factors.const_ref();
        af::shared<double> result(uffs.size(), af::init_functor_null<double>());
        double x_sq = d_star_sq / 4;
        for(std::size_t i=0;i<uffs.size();i++) {
          if (!uffs[i]) {
            throw std::runtime_error(
              "gaussian not defined for scattering_type \""
              + type_given_unique_index(i)
              + "\".");
          }
          result[i] = uffs[i]->at_x_sq(x_sq);
        }
        return result;
      }
  };

  typedef generic_scattering_type_registry<eltbx::xray_scattering::gaussian>
          scattering_type_registry;


  template<>
  struct form_factor_table_policy<eltbx::xray_scattering::gaussian>
  {
    typedef eltbx::xray_scattering::gaussian form_factor_t;

    template<class TypeIndexPairsConstIter>
    static void
    assign_from(std::string const &table,
                af::ref<boost::optional<form_factor_t> > const &unique_ffs,
                TypeIndexPairsConstIter first_type_index_pair,
                TypeIndexPairsConstIter last_type_index_pair)
    {
      CCTBX_ASSERT(table == "IT1992" || table == "WK1995");
      if (table == "IT1992") {
        for(TypeIndexPairsConstIter
            pair=first_type_index_pair;
            pair!=last_type_index_pair;
            pair++) {
          std::size_t ui = pair->second;
          if (unique_ffs[ui]) continue;
          unique_ffs[ui] = eltbx::xray_scattering::it1992(pair->first,
                                                          true).fetch();
        }
      }
      else {
        for(TypeIndexPairsConstIter
            pair=first_type_index_pair;
            pair!=last_type_index_pair;
            pair++) {
          std::size_t ui = pair->second;
          if (unique_ffs[ui]) continue;
          unique_ffs[ui] = eltbx::xray_scattering::wk1995(pair->first,
                                                          true).fetch();
        }
      }
    }

  };

}} // namespace cctbx::xray

#endif // CCTBX_XRAY_SCATTERING_TYPE_REGISTRY_H
