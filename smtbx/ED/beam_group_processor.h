#pragma once
#include <cctbx/xray/thickness.h>

#include <smtbx/refinement/least_squares_fc.h>
#include <smtbx/ED/utils.h>
#include <smtbx/ED/ed_data.h>
#include <smtbx/refinement/constraints/reparametrisation.h>
#include <scitbx/array_family/selections.h>
#include <smtbx/ED/dyn_calculator.h>

namespace smtbx {  namespace ED
{
  using namespace refinement::least_squares;
  template <typename FloatType>
  class beam_group_processor {
  public:
    ED_UTIL_TYPEDEFS;
    typedef boost::shared_ptr<a_dyn_calculator<FloatType> > dyn_calculator_ptr_t;
    beam_group_processor(
      const dyn_calculator_factory<FloatType>& dc_f,
      const BeamGroup<FloatType>& beam_group,
      const cmat_t& Ugs,
      const cart_t& K,
      const cctbx::xray::thickness<FloatType>& thickness,
      // kinematic dFc_dP
      const af::shared<cmat_t>& Ds_kin,
      bool calc_grad)
      : beam_group(beam_group),
      Ugs(Ugs),
      thickness(thickness),
      Ds_kin(Ds_kin),
      calc_grad(calc_grad)
    {
      strong_indices = af::select(beam_group.indices.const_ref(),
        beam_group.strong_beams.const_ref());
      dyn_calculator = dc_f.make(strong_indices, K, thickness.value);
    }

    virtual ~beam_group_processor() {}

    void operator()() {
      process(beam_group.RMf, beam_group.geometry.get_normal());
    }

    void process(const mat3_t& RM, const cart_t& N) {
      try {
        size_t beam_n = beam_group.strong_measured_beams.size();
        if (calc_grad) {
          CIs = dyn_calculator->reset(Ugs, RM, N).calc_amps_ext(Ds_kin,
            thickness.grad,
            D_dyn,
            beam_n);
        }
        else {
          CIs = dyn_calculator->reset(Ugs, RM, N)
            .calc_amps(beam_n);
        }
      }
      catch (smtbx::error const& e) {
        exception_.reset(new smtbx::error(e));
      }
      catch (std::exception const& e) {
        exception_.reset(new smtbx::error(e.what()));
      }
    }

    void process_1(size_t idx, const mat3_t& RM, const cart_t& N) {
      try {
        if (calc_grad) {
          CIs.resize(1);
          CIs[0] = dyn_calculator->reset(Ugs, RM, N).calc_amps_ext_1(Ds_kin,
            thickness.grad,
            D_dyn,
            idx);
        }
        else {
          CIs.resize(1);
          CIs[0] = dyn_calculator->reset(Ugs, RM, N).calc_amps_1(idx);
        }
      }
      catch (smtbx::error const& e) {
        exception_.reset(new smtbx::error(e));
      }
      catch (std::exception const& e) {
        exception_.reset(new smtbx::error(e.what()));
      }
    }

    const cart_t& getK() const {
      return dyn_calculator->K;
    }
  private:
    dyn_calculator_ptr_t dyn_calculator;
    af::shared<miller::index<> > strong_indices;
  public:
    const BeamGroup<FloatType>& beam_group;
    cmat_t Ugs;
    const cctbx::xray::thickness<FloatType>& thickness;
    af::shared<cmat_t> Ds_kin;
    const bool calc_grad;
    boost::scoped_ptr<smtbx::error> exception_;
    // results
    af::shared<complex_t> CIs;
    mat_t D_dyn;
  };

}}
