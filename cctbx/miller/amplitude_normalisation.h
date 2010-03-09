#ifndef CCTBX_MILLER_NORMALISED_ARRAY_H
#define CCTBX_MILLER_NORMALISED_ARRAY_H

#include <cctbx/uctbx.h>
#include <cctbx/sgtbx/space_group.h>
#include <cctbx/miller.h>
#include <cctbx/eltbx/xray_scattering/gaussian.h>
#include <cctbx/import_scitbx_af.h>

namespace cctbx { namespace miller {

  /// Amplitudes normalisation transforming so-called F's into E's, and vice-versa.
  template <typename FloatType>
  class amplitude_normalisation
  {
  public:
    typedef eltbx::xray_scattering::gaussian form_factor_t;
    typedef FloatType float_type;

    amplitude_normalisation(af::const_ref<form_factor_t> const &form_factors,
                            af::const_ref<float_type> const &multiplicities,
                            float_type wilson_intensity_scale_factor,
                            float_type wilson_b,
                            uctbx::unit_cell const &unit_cell,
                            sgtbx::space_group const &space_group,
                            af::const_ref<index<> > const &indices)
    {
      CCTBX_ASSERT(form_factors.size() == multiplicities.size())
                  (form_factors.size())(multiplicities.size());
      std::size_t n = indices.size();
      normalisations.reserve(n);
      for (std::size_t i=0; i<n; ++i) {
        index<> const &h = indices[i];
        float_type stol_sq = unit_cell.stol_sq(h);
        float_type sum_form_factor_sq = 0;
        for (std::size_t j = 0; j != form_factors.size(); j++)
        {
          float_type n_atoms = multiplicities[j] * space_group.order_z();
          float_type form_factor = form_factors[j].at_stol_sq(stol_sq);
          sum_form_factor_sq += n_atoms*form_factor*form_factor;
        }
        float_type normalisation_at_h = wilson_intensity_scale_factor
                                          * std::exp(-2 * wilson_b * stol_sq)
                                          * space_group.epsilon(h)
                                          * space_group.n_ltr()
                                          * sum_form_factor_sq;
        normalisations.push_back(std::sqrt(normalisation_at_h));
      }
    }

    /// E(h) = F(h) / N(h) where the N's are the elements of that array.
    af::shared<float_type> normalisations;
  };

}} // cctbx::miller

#endif
