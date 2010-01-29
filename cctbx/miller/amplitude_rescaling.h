#ifndef CCTBX_MILLER_NORMALISED_ARRAY_H
#define CCTBX_MILLER_NORMALISED_ARRAY_H

#include <cctbx/uctbx.h>
#include <cctbx/sgtbx/space_group.h>
#include <cctbx/miller.h>
#include <cctbx/eltbx/xray_scattering/gaussian.h>
#include <cctbx/import_scitbx_af.h>

namespace cctbx { namespace miller {

  //! A class for normalising a set of amplitudes, or the reverse operation.
  template <typename FloatType>
  class amplitude_rescaling
  {
  public:
    typedef eltbx::xray_scattering::gaussian form_factor_t;

    enum rescaling_kind { normalised, denormalised };

    amplitude_rescaling(
      int kind, // {normalised, denormalised}
      af::const_ref<form_factor_t> const &form_factors_,
      af::const_ref<FloatType> const &multiplicities_,
      FloatType wilson_intensity_scale_factor_,
      FloatType wilson_b_,
      uctbx::unit_cell const &unit_cell_,
      sgtbx::space_group const &space_group_,
      af::shared<index<> > const &indices,
      af::shared<FloatType> const &data)
    : form_factors(form_factors_),
      multiplicities(multiplicities_),
      wilson_intensity_scale_factor(wilson_intensity_scale_factor_),
      wilson_b(wilson_b_),
      unit_cell(unit_cell_),
      space_group(space_group_),
      sum_e_sq_minus_1(0),
      n_e_greater_than_2(0)
    {
      CCTBX_ASSERT(kind == normalised || kind == denormalised);
      CCTBX_ASSERT(form_factors.size() == multiplicities.size())
                  (form_factors.size())(multiplicities.size());

      for(std::size_t i=0;i<indices.size();i++) {
        FloatType f_sq, e_sq;
        index<> const &h = indices[i];
        if (kind == normalised) {
          f_sq = std::pow(data[i], 2);
          e_sq = f_sq / normalising_denominator(h);
          result.push_back(std::sqrt(e_sq));
        }
        else if (kind == denormalised) {
          e_sq = std::pow(data[i], 2);
          f_sq = e_sq * normalising_denominator(h);
          result.push_back(std::sqrt(f_sq));
        }
        sum_e_sq_minus_1 += std::abs(e_sq - 1);
        if (e_sq > 4.0) n_e_greater_than_2 += 1;
      }
    }

    FloatType sum_e_sq_minus_1;
    int n_e_greater_than_2;
    af::shared<FloatType> result;

  private:
    // The denominator in the equation to normalise a given reflection
    FloatType normalising_denominator(index<> const &h) {
      FloatType stol_sq = unit_cell.stol_sq(h);
      int epsilon = space_group.epsilon(h);
      FloatType sum_form_factor_sq = 0;
      for (std::size_t j = 0; j != form_factors.size(); j++)
      {
        FloatType n_atoms = multiplicities[j] * space_group.order_z();
        FloatType form_factor = form_factors[j].at_stol_sq(stol_sq);
        sum_form_factor_sq += n_atoms*form_factor*form_factor;
      }
      return wilson_intensity_scale_factor
           * std::exp(-2 * wilson_b * stol_sq)
           * space_group.epsilon(h)
           * space_group.n_ltr()
           * sum_form_factor_sq;
    }

    af::const_ref<form_factor_t> const &form_factors;
    af::const_ref<FloatType> const &multiplicities;
    FloatType wilson_intensity_scale_factor;
    FloatType wilson_b;
    uctbx::unit_cell const &unit_cell;
    sgtbx::space_group const &space_group;

  };

}} // cctbx::miller

#endif
