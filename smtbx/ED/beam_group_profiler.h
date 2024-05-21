#pragma once

#include <cctbx/xray/thickness.h>
#include <smtbx/ED/utils.h>
#include <smtbx/ED/ed_data.h>
#include <smtbx/ED/dyn_calculator.h>
#include <smtbx/ED/n_beam.h>
#include <smtbx/refinement/least_squares_fc.h>
#include <cctbx/xray/fc_correction.h>
#include <boost/scoped_ptr.hpp>

namespace smtbx { namespace ED {
  using namespace smtbx::refinement::least_squares;

  template <typename FloatType>
  struct beam_group_profiler {
    ED_UTIL_TYPEDEFS;
    typedef f_calc_function_base<FloatType> f_calc_function_base_t;
    typedef boost::shared_ptr< f_calc_function_base_t> f_calc_function_base_ptr_t;
    typedef af::shared<const BeamInfo<FloatType>*> beam_at;
    typedef std::pair<int, af::shared<const BeamInfo<FloatType>*> > beam_me;
    typedef boost::shared_ptr<lookup_t> lookup_ptr_t;
    typedef cctbx::xray::fc_correction<FloatType> fc_correction_t;
    typedef boost::shared_ptr< fc_correction_t> fc_correction_ptr_t;

    beam_group_profiler(const BeamGroup<FloatType>& beam_group,
      f_calc_function_base_t& f_calc_function,
      const fc_correction_t& fc_cr,
      const sgtbx::space_group& space_group,
      bool anomalous_flag,
      const cctbx::xray::thickness<FloatType>& thickness,
      const RefinementParams<FloatType>& params)
      : beam_group(beam_group),
      f_calc_function(f_calc_function),
      space_group(space_group),
      fc_cr(fc_cr),
      params(params),
      Kl_(params.getKl()),
      Fc2Ug(params.getFc2Ug()),
      mat_type(params.getMatrixType()),
      thickness(thickness.value),
      thread_n(params.getThreadN())
    {
      if (mat_type >= 100) {
        mat_type -= 100;
        use_n_beam = true;
      }
      else {
        use_n_beam = false;
      }
      K = beam_group.geometry->Kl_as_K(Kl_);
      // build lookups for each beam_group + collect all indices and they diffs
      af::shared<miller::index<> > all_indices;
      // treat equivalents independently inside the beam_groups
      sgtbx::space_group P1("P 1");
      for (size_t hi = 0; hi < beam_group.indices.size(); hi++) {
        all_indices.push_back(beam_group.indices[hi]);
        for (size_t hj = hi + 1; hj < beam_group.indices.size(); hj++) {
          all_indices.push_back(beam_group.indices[hi] - beam_group.indices[hj]);
          all_indices.push_back(beam_group.indices[hj] - beam_group.indices[hi]);
        }
      }
      mi_lookup = lookup_t(
        all_indices.const_ref(),
        P1,
        anomalous_flag);
      indices = mi_lookup.get_unique();
      mi_lookup = lookup_t(
        indices.const_ref(),
        P1,
        anomalous_flag);
      strong_indices = af::select(beam_group.indices.const_ref(),
        beam_group.strong_beams.const_ref());
      build_();
      if (!use_n_beam) {
        utils<FloatType>::build_Ug_matrix(
          A, Fcs_k,
          mi_lookup,
          strong_indices);
      }
    }

    complex_t calc_one_h(miller::index<> const& h) const {
      f_calc_function.compute(h, boost::none, 0, false);
      complex_t fc = f_calc_function.get_f_calc();
      FloatType fc_k = fc_cr.compute(h, f_calc_function.get_observable(), false);
      if (fc_k != 1) {
        fc *= std::sqrt(fc_k);
      }
      return fc;
    }

    struct process_beam_group_profile {
      process_beam_group_profile(const beam_group_profiler& parent,
        const BeamGroup<FloatType>& beam_group)
        : parent(parent),
        beam_group(beam_group)
      {}

      // N-Beam specific
      af::shared<FloatType> process(const miller::index<> &h, const af::shared<FloatType> &angles) {
        SMTBX_ASSERT(parent.use_n_beam && parent.params.getBeamN() > 2);
        try {
          af::shared<FloatType> rv(angles.size());
          bool floating = parent.params.isNBeamFloating();
          dyn_calculator_n_beam<FloatType> n_beam_dc(parent.params.getBeamN(),
            parent.mat_type,
            beam_group, parent.K, parent.thickness,
            parent.params.useNBeamSg(), parent.params.getNBeamWght());
          if (!floating) {
            FloatType da = beam_group.get_diffraction_angle(h, parent.K);
            n_beam_dc.init(h, da, parent.Fcs_k, parent.mi_lookup);
          }
          for (size_t ai = 0; ai < angles.size(); ai++) {
            std::pair<mat3_t, cart_t> r = beam_group.compute_RMf_N(angles[ai]);
            if (floating) {
              n_beam_dc.init(h, r.first, parent.Fcs_k, parent.mi_lookup);
            }
            rv[ai] = std::norm(n_beam_dc.calc_amp(r));
          }
          return rv;
        }
        catch (smtbx::error const& e) {
          exception_.reset(new smtbx::error(e));
        }
        catch (std::exception const& e) {
          exception_.reset(new smtbx::error(e.what()));
        }
      }

      af::shared<FloatType> process(const mat3_t& RMf, const cart_t& N) {
        af::shared<FloatType> Is(beam_group.strong_measured_beams.size());
        try {
          if (parent.use_n_beam) { // N-beam
            /* a LOT of overhead here!!! may need to change the logic to speed up
             as the matrix rebuilt for each angle
             */
            if (parent.params.getBeamN() > 2) {
              bool floating = parent.params.isNBeamFloating();
              for (size_t i = 0; i < beam_group.strong_measured_beams.size(); i++) {
                size_t beam_idx = beam_group.strong_measured_beams[i];
                miller::index<> h = beam_group.indices[beam_idx];
                dyn_calculator_n_beam<FloatType> n_beam_dc(parent.params.getBeamN(),
                  parent.mat_type,
                  beam_group, parent.K, parent.thickness,
                  parent.params.useNBeamSg(), parent.params.getNBeamWght());
                FloatType da = beam_group.get_diffraction_angle(h, parent.K);
                std::pair<mat3_t, cart_t> da_r = beam_group.compute_RMf_N(da);
                if (floating) {
                  n_beam_dc.init(h, RMf, parent.Fcs_k, parent.mi_lookup);
                }
                else {
                  n_beam_dc.init(h, da_r.first, parent.Fcs_k, parent.mi_lookup);
                }
                Is[i] = std::norm(
                  n_beam_dc.calc_amp(std::make_pair(RMf, da_r.second)));
              }
            }
            else {
              for (size_t i = 0; i < beam_group.strong_measured_beams.size(); i++) {
                size_t beam_idx = beam_group.strong_measured_beams[i];
                miller::index<> h = beam_group.indices[beam_idx];
                int ii = parent.mi_lookup.find_hkl(h);
                complex_t Fc = ii != -1 ? parent.Fcs_k[ii] : 0;
                FloatType da = beam_group.get_diffraction_angle(h, parent.K);
                std::pair<mat3_t, cart_t> da_r = beam_group.compute_RMf_N(da);
                Is[i] = std::norm(utils<FloatType>::calc_amp_2beam(
                  h, Fc,
                  parent.thickness,
                  parent.K, RMf, da_r.second));
              }
            }
          }
          else {
            boost::shared_ptr<a_dyn_calculator<FloatType> > dc =
              dyn_calculator_factory<FloatType>(parent.mat_type)
              .make(parent.strong_indices, parent.K, parent.thickness);
            af::shared<complex_t> amps =
              dc->reset(parent.A, RMf, beam_group.geometry->get_normal()).calc_amps(
                beam_group.strong_measured_beams.size());
            for (size_t i = 0; i < amps.size(); i++) {
              Is[i] = std::norm(amps[i]);
            }
          }
        }
        catch (smtbx::error const& e) {
          exception_.reset(new smtbx::error(e));
        }
        catch (std::exception const& e) {
          exception_.reset(new smtbx::error(e.what()));
        }
        return Is;
      }
      
      FloatType process_incident(const mat3_t& RMf, const cart_t& N) {
        try {
          if (parent.use_n_beam) { // N-beam
            FloatType I_sum = 0;
            if (parent.params.getBeamN() > 2) {
              bool floating = parent.params.isNBeamFloating();
              for (size_t i = 0; i < beam_group.strong_measured_beams.size(); i++) {
                size_t beam_idx = beam_group.strong_measured_beams[i];
                miller::index<> h = beam_group.indices[beam_idx];
                dyn_calculator_n_beam<FloatType> n_beam_dc(parent.params.getBeamN(),
                  parent.mat_type,
                  beam_group, parent.K, parent.thickness,
                  parent.params.useNBeamSg(), parent.params.getNBeamWght());
                if (floating) {
                  n_beam_dc.init(h, RMf, parent.Fcs_k, parent.mi_lookup);
                }
                else {
                  FloatType da = beam_group.get_diffraction_angle(h, parent.K);
                  n_beam_dc.init(h, da, parent.Fcs_k, parent.mi_lookup);
                }
                I_sum += std::norm(
                  n_beam_dc.calc_amp(std::make_pair(RMf, N), 0));
              }
            }
            else {
              for (size_t i = 0; i < beam_group.strong_measured_beams.size(); i++) {
                size_t beam_idx = beam_group.strong_measured_beams[i];
                miller::index<> h = beam_group.indices[beam_idx];
                int ii = parent.mi_lookup.find_hkl(h);
                complex_t Fc = ii != -1 ? parent.Fcs_k[ii] : 0;
                I_sum += std::norm(utils<FloatType>::calc_amp_2beam(
                  h, Fc,
                  parent.thickness,
                  parent.K, RMf, N, 0));
              }
            }
            return I_sum / beam_group.strong_measured_beams.size();
          }
          boost::shared_ptr<a_dyn_calculator<FloatType> > dc =
            dyn_calculator_factory<FloatType>(parent.mat_type)
            .make(parent.strong_indices, parent.K, parent.thickness);
          return std::norm(dc->reset(parent.A, RMf, N).calc_amps_1(0));
        }
        catch (smtbx::error const& e) {
          exception_.reset(new smtbx::error(e));
        }
        catch (std::exception const& e) {
          exception_.reset(new smtbx::error(e.what()));
        }
        // suppress warning
        return FloatType();
      }
      boost::shared_ptr< dyn_calculator_n_beam<FloatType> > n_beam_dc;
      const beam_group_profiler& parent;
      const BeamGroup<FloatType>& beam_group;
      FloatType thickness;
      boost::scoped_ptr<smtbx::error> exception_;
    };

    af::shared<PeakProfilePoint<FloatType> > build_profile(
      const af::shared<FloatType>& angles, bool inc_incident)
    {
      af::shared<PeakProfilePoint<FloatType> > rv;
      for (size_t ai = 0; ai < angles.size(); ai++) {
        FloatType ang = angles[ai];
        std::pair<mat3_t, cart_t> r = beam_group.compute_RMf_N(ang);
        process_beam_group_profile fp(*this, beam_group);
        af::shared<FloatType> Is = fp.process(r.first, r.second);
        if (fp.exception_) {
          throw* fp.exception_.get();
        }
        if (inc_incident) {
          process_beam_group_profile ibp(*this, beam_group);
          FloatType I = ibp.process_incident(r.first, r.second);
          rv.push_back(PeakProfilePoint<FloatType>(I, 0.0, ang, Kl_));
        }
        for (size_t i = 0; i < beam_group.strong_measured_beams.size(); i++) {
          size_t beam_idx = beam_group.strong_measured_beams[i];
          miller::index<> h = beam_group.indices[beam_idx];
          cart_t g = r.first * cart_t(h[0], h[1], h[2]);
          FloatType Sg = utils<FloatType>::calc_Sg(g, K);
          rv.push_back(PeakProfilePoint<FloatType>(Is[i], Sg, ang, (K + g).length()));
        }
      }
      return rv;
    }

    af::shared<PeakProfilePoint<FloatType> > build_incident_profile(
      const af::shared<FloatType>& angles)
    {
      af::shared<PeakProfilePoint<FloatType> > rv(
        af::reserve(angles.size()));
      for (size_t ai = 0; ai < angles.size(); ai++) {
        FloatType ang = angles[ai];
        std::pair<mat3_t, cart_t> r = beam_group.compute_RMf_N(ang);
        process_beam_group_profile fp(*this, beam_group);
        FloatType I = fp.process_incident(r.first, r.second);
        if (fp.exception_) {
          throw* fp.exception_.get();
        }
        rv.push_back(PeakProfilePoint<FloatType>(I, 0.0, ang, Kl_));
      }
      return rv;
    }

    void build_() {
      Fcs_k.resize(indices.size());
      for (size_t ih = 0; ih < indices.size(); ih++) {
        Fcs_k[ih] = calc_one_h(indices[ih]) * Fc2Ug;
      }
    }

    const BeamGroup<FloatType>& beam_group;
    f_calc_function_base_t& f_calc_function;
    const cctbx::xray::fc_correction<FloatType>& fc_cr;
    const sgtbx::space_group& space_group;
    af::shared<miller::index<> > indices;
    RefinementParams<FloatType> params;
    FloatType Fc2Ug, Kl_;
    cart_t K;
    int mat_type;
    size_t beam_n;
    FloatType thickness;
    af::shared<complex_t> Fcs_k;
    lookup_t mi_lookup;
    int thread_n;
    bool use_n_beam;
    af::shared<miller::index<> > strong_indices;
    cmat_t A;
  };
}} // end of smtbx::ED
