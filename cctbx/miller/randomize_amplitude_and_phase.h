#ifndef CCTBX_RANDOMIZE_AMPLITUDE_AND_PHASE_H
#define CCTBX_RANDOMIZE_AMPLITUDE_AND_PHASE_H

#include <cctbx/import_scitbx_af.h>
#include <cctbx/error.h>
#include <scitbx/random.h>

namespace cctbx { namespace miller {

  //Randomize complex number with prescribed target error in amplitude and phase
  template <typename ComplexType, typename FloatType>
  af::shared<ComplexType> randomize_amplitude_and_phase(
    af::const_ref<ComplexType> const& data,
    af::const_ref<bool> const& selection,
    FloatType const& amplitude_error,
    FloatType const& phase_error_deg,
    int random_seed)
  {
    CCTBX_ASSERT(amplitude_error>=0 && amplitude_error<=1);
    CCTBX_ASSERT(phase_error_deg>=0 && phase_error_deg<=90);
    CCTBX_ASSERT(data.size() == selection.size());
    scitbx::random::mersenne_twister mt(random_seed);
    af::shared<ComplexType>
      result(data.size(), af::init_functor_null<ComplexType>());
    for (std::size_t i=0; i<data.size(); i++) {
      if(selection[i]) {
        ComplexType d = data[i];
        FloatType a = std::abs(d);
        FloatType p = std::arg(d);
        FloatType sc = 1;
        if(mt.random_double() <= 0.5) sc=-1;
        //sc *= (2 * mt.random_double()); //Why only works up to threshold=0.5?
        a += a*amplitude_error*sc;
        p += p*(phase_error_deg/90.)*sc;
        result[i] = ComplexType(a*std::cos(p), a*std::sin(p));
      }
      else {
        result[i] = data[i];
      }
    }
    return result;
  };

}} // namespace cctbx::miller

#endif // CCTBX_RANDOMIZE_AMPLITUDE_AND_PHASE_H
