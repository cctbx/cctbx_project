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
  struct frame_profiler {
    ED_UTIL_TYPEDEFS;
    typedef f_calc_function_base<FloatType> f_calc_function_base_t;
    typedef boost::shared_ptr< f_calc_function_base_t> f_calc_function_base_ptr_t;
    typedef af::shared<const BeamInfo<FloatType>*> beam_at;
    typedef std::pair<int, af::shared<const BeamInfo<FloatType>*> > beam_me;
    typedef boost::shared_ptr<lookup_t> lookup_ptr_t;
    typedef cctbx::xray::fc_correction<FloatType> fc_correction_t;
    typedef boost::shared_ptr< fc_correction_t> fc_correction_ptr_t;

    frame_profiler(const FrameInfo<FloatType>& frame,
      f_calc_function_base_t& f_calc_function,
      const fc_correction_t& fc_cr,
      const sgtbx::space_group& space_group,
      bool anomalous_flag,
      const cctbx::xray::thickness<FloatType>& thickness,
      const RefinementParams<FloatType>& params)
      : frame(frame),
      f_calc_function(f_calc_function),
      space_group(space_group),
      fc_cr(fc_cr),
      params(params),
      Kl(params.getKl()),
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
      // build lookups for each frame + collect all indices and they diffs
      af::shared<miller::index<> > all_indices;
      // treat equivalents independently inside the frames
      sgtbx::space_group P1("P 1");
      for (size_t hi = 0; hi < frame.indices.size(); hi++) {
        all_indices.push_back(frame.indices[hi]);
        for (size_t hj = hi + 1; hj < frame.indices.size(); hj++) {
          all_indices.push_back(frame.indices[hi] - frame.indices[hj]);
          all_indices.push_back(frame.indices[hj] - frame.indices[hi]);
        }
      }
      mi_lookup = lookup_t(
        all_indices.const_ref(),
        space_group,
        anomalous_flag);
      indices = mi_lookup.get_unique();
      mi_lookup = lookup_t(
        indices.const_ref(),
        space_group,
        anomalous_flag);
      strong_indices = af::select(frame.indices.const_ref(),
        frame.strong_beams.const_ref());
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

    struct process_frame_profile {
      process_frame_profile(const frame_profiler& parent,
        const FrameInfo<FloatType>& frame,
        af::shared<FloatType>& Is)
        : parent(parent),
        frame(frame),
        Is(Is)
      {
      }

      void process(const mat3_t& RMf, const cart_t& N) {
        try {
          const cart_t K = cart_t(0, 0, -parent.Kl);
          if (parent.use_n_beam) { // N-beam
            if (parent.params.getBeamN() > 2) {
              for (size_t i = 0; i < frame.strong_measured_beams.size(); i++) {
                size_t beam_idx = frame.strong_measured_beams[i];
                miller::index<> h = frame.indices[beam_idx];
                dyn_calculator_n_beam<FloatType> n_beam_dc(parent.params.getBeamN(),
                  parent.mat_type,
                  frame, K, parent.thickness);
                n_beam_dc.init(h, RMf, parent.Fcs_k, parent.mi_lookup);
                Is[i] = std::norm(
                  n_beam_dc.calc_amp(std::make_pair(RMf, N)));
              }
            }
            else {
              for (size_t i = 0; i < frame.strong_measured_beams.size(); i++) {
                size_t beam_idx = frame.strong_measured_beams[i];
                miller::index<> h = frame.indices[beam_idx];
                int ii = parent.mi_lookup.find_hkl(h);
                complex_t Fc = ii != -1 ? parent.Fcs_k[ii] : 0;
                Is[i] = std::norm(utils<FloatType>::calc_amp_2beam(
                  h, Fc,
                  parent.thickness,
                  K, RMf, N));
              }
            }
            return;
          }
          boost::shared_ptr<a_dyn_calculator<FloatType> > dc =
            dyn_calculator_factory<FloatType>(parent.mat_type)
              .make(parent.strong_indices, K, parent.thickness);
          af::shared<complex_t> amps = 
            dc->reset(parent.A, RMf, N).calc_amps(
              frame.strong_measured_beams.size());
          for (size_t i = 0; i < amps.size(); i++) {
            Is[i] = std::norm(amps[i]);
          }
        }
        catch (smtbx::error const& e) {
          exception_.reset(new smtbx::error(e));
        }
        catch (std::exception const& e) {
          exception_.reset(new smtbx::error(e.what()));
        }
      }
      
      FloatType process_incident(const mat3_t& RMf, const cart_t& N) {
        try {
          const cart_t K = cart_t(0, 0, -parent.Kl);
          if (parent.use_n_beam) { // N-beam
            FloatType I_sum = 0;
            if (parent.params.getBeamN() > 2) {
              for (size_t i = 0; i < frame.strong_measured_beams.size(); i++) {
                size_t beam_idx = frame.strong_measured_beams[i];
                miller::index<> h = frame.indices[beam_idx];
                dyn_calculator_n_beam<FloatType> n_beam_dc(parent.params.getBeamN(),
                  parent.mat_type,
                  frame, K, parent.thickness);
                n_beam_dc.init(h, RMf, parent.Fcs_k, parent.mi_lookup);
                I_sum += std::norm(
                  n_beam_dc.calc_amp(std::make_pair(RMf, N), 0));
              }
            }
            else {
              for (size_t i = 0; i < frame.strong_measured_beams.size(); i++) {
                size_t beam_idx = frame.strong_measured_beams[i];
                miller::index<> h = frame.indices[beam_idx];
                int ii = parent.mi_lookup.find_hkl(h);
                complex_t Fc = ii != -1 ? parent.Fcs_k[ii] : 0;
                I_sum += std::norm(utils<FloatType>::calc_amp_2beam(
                  h, Fc,
                  parent.thickness,
                  K, RMf, N, 0));
              }
            }
            return I_sum / frame.strong_measured_beams.size();
          }
          boost::shared_ptr<a_dyn_calculator<FloatType> > dc =
            dyn_calculator_factory<FloatType>(parent.mat_type)
            .make(parent.strong_indices, K, parent.thickness);
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
      const frame_profiler& parent;
      const FrameInfo<FloatType>& frame;
      FloatType thickness;
      af::shared<FloatType>& Is;
      boost::scoped_ptr<smtbx::error> exception_;
    };

    af::shared<PeakProfilePoint<FloatType> > build_profile(
      const af::shared<FloatType>& angles, bool inc_incident)
    {
      af::shared<FloatType> Is_(frame.beams.size());
      af::shared<PeakProfilePoint<FloatType> > rv;
      const cart_t K = cart_t(0, 0, -Kl);
      for (size_t ai = 0; ai < angles.size(); ai++) {
        FloatType ang = angles[ai];
        std::pair<mat3_t, cart_t> r = frame.compute_RMf_N(ang);
        std::fill(Is_.begin(), Is_.end(), 0);
        process_frame_profile fp(*this, frame, Is_);
        fp.process(r.first, r.second);
        if (fp.exception_) {
          throw* fp.exception_.get();
        }
        if (inc_incident) {
          process_frame_profile ibp(*this, frame, Is_);
          FloatType I = ibp.process_incident(r.first, r.second);
          rv.push_back(PeakProfilePoint<FloatType>(I, 0.0, ang, Kl));
        }
        for (size_t i = 0; i < frame.strong_measured_beams.size(); i++) {
          size_t beam_idx = frame.strong_measured_beams[i];
          miller::index<> h = frame.indices[beam_idx];
          cart_t g = r.first * cart_t(h[0], h[1], h[2]);
          FloatType Sg = utils<FloatType>::calc_Sg(g, K);
          rv.push_back(PeakProfilePoint<FloatType>(Is_[i], Sg, ang, (K + g).length()));
        }
      }
      return rv;
    }

    af::shared<PeakProfilePoint<FloatType> > build_incident_profile(
      const af::shared<FloatType>& angles)
    {
      af::shared<FloatType> Is_;
      af::shared<PeakProfilePoint<FloatType> > rv(
        af::reserve(angles.size()));
      const cart_t K = cart_t(0, 0, -Kl);
      for (size_t ai = 0; ai < angles.size(); ai++) {
        FloatType ang = angles[ai];
        std::pair<mat3_t, cart_t> r = frame.compute_RMf_N(ang);
        process_frame_profile fp(*this, frame, Is_);
        FloatType I = fp.process_incident(r.first, r.second);
        if (fp.exception_) {
          throw* fp.exception_.get();
        }
        rv.push_back(PeakProfilePoint<FloatType>(I, 0.0, ang, Kl));
      }
      return rv;
    }

    void build_() {
      Fcs_k.resize(indices.size());
      for (size_t ih = 0; ih < indices.size(); ih++) {
        Fcs_k[ih] = calc_one_h(indices[ih]) * Fc2Ug;
      }
    }

    const FrameInfo<FloatType>& frame;
    f_calc_function_base_t& f_calc_function;
    const cctbx::xray::fc_correction<FloatType>& fc_cr;
    const sgtbx::space_group& space_group;
    af::shared<miller::index<> > indices;
    RefinementParams<FloatType> params;
    FloatType Kl, Fc2Ug;
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
