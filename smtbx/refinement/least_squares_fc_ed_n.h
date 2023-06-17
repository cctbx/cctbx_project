#ifndef SMTBX_REFINEMENT_LEAST_SQUARES_FC_ED_N_H
#define SMTBX_REFINEMENT_LEAST_SQUARES_FC_ED_N_H

#include <cctbx/xray/thickness.h>

#include <smtbx/refinement/least_squares_fc.h>
#include <smtbx/ED/ed_data.h>
#include <smtbx/refinement/constraints/reparametrisation.h>
#include <scitbx/array_family/selections.h>

namespace smtbx {  namespace refinement {
namespace least_squares {
  using namespace smtbx::ED;
  using namespace refinement::constraints;

  template <typename FloatType>
  struct ed_n_shared_data {
    typedef std::complex<FloatType> complex_t;
    typedef f_calc_function_base<FloatType> f_calc_function_base_t;
    typedef boost::shared_ptr< f_calc_function_base_t> f_calc_function_base_ptr_t;
    typedef builder_base<FloatType> data_t;
    typedef af::shared<const BeamInfo<FloatType>*> beam_at;
    typedef std::pair<int, af::shared<const BeamInfo<FloatType>*> > beam_me;
    typedef miller::lookup_utils::lookup_tensor<FloatType> lookup_t;
    typedef boost::shared_ptr<lookup_t> lookup_ptr_t;
    typedef cctbx::xray::fc_correction<FloatType> fc_correction_t;
    typedef boost::shared_ptr< fc_correction_t> fc_correction_ptr_t;

    ed_n_shared_data(reparametrisation const& reparamn,
      f_calc_function_base_t& f_calc_function,
      fc_correction_t const& fc_cr,
      sgtbx::space_group const& space_group,
      bool anomalous_flag,
      af::shared<FrameInfo<FloatType> > frames,
      cctbx::xray::thickness<FloatType> const& thickness,
      // Kl, Fc2Ug, epsilon
      af::shared<FloatType> const& params,
      bool compute_grad,
      bool do_build = true)
      : reparamn(reparamn),
      f_calc_function(f_calc_function),
      space_group(space_group),
      fc_cr(fc_cr),
      Kvac(params[0]),
      Kl(params[1]),
      Fc2Ug(params[2]),
      eps(params[3]),
      mat_type(static_cast<int>(params[4])),
      frames(frames),
      thickness(thickness),
      compute_grad(compute_grad),
      thread_n(static_cast<int>(params[5])),
      use_diff_angle(params[6] > 0)
    {
      // build lookups for each frame + collect all indices and they diffs
      af::shared<miller::index<> > all_indices;
      size_t offset = 0;
      // treat equivalents independently inside the frames
      sgtbx::space_group P1("P 1");
      for (size_t i=0; i < frames.size(); i++) {
        FrameInfo<FloatType>& frame = frames[i];
        frames_map.insert(std::make_pair(frame.id, &frame));
        frame.offset = offset;
        offset += frame.strong_measured_beams.size();
        for (size_t hi = 0; hi < frame.indices.size(); hi++) {
          all_indices.push_back(frame.indices[hi]);
          for (size_t hj = hi+1; hj < frame.indices.size(); hj++) {
            all_indices.push_back(frame.indices[hi] - frame.indices[hj]);
          }
        }
        lookup_ptr_t mi_l(new lookup_t(
          af::select(frame.indices.const_ref(),
            frame.strong_measured_beams.const_ref()).const_ref(),
          P1,
          false));
        frame_lookups.insert(std::make_pair(frames[i].id, mi_l));
      }
      beam_n = offset;
      // a tricky way of getting unique only...
      mi_lookup = lookup_t(
        all_indices.const_ref(),
        space_group,
        anomalous_flag);
      indices = mi_lookup.get_unique();
      mi_lookup = lookup_t(
        indices.const_ref(),
        space_group,
        anomalous_flag);
      if (do_build) {
        build();
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

    /* could be used for both - values and derivatives calculation for the given
    set of beams of a frame. The result is written to Fcs starting at the given
    offset
    */
    struct process_frame {
      process_frame(ed_n_shared_data const& parent,
        FrameInfo<FloatType> &frame,
        // source kinematic Fcs, Fc, Fc_+e, Fc_-e
        const af::shared<complex_t> Fcs_k,
        af::shared<FloatType>& Is,
        af::shared<complex_t>& CIs, bool use_offset,
        bool use_diff_angle)
        : parent(parent),
        frame(frame),
        thickness(parent.thickness.value),
        Fcs_k(Fcs_k),
        Is(Is),
        CIs(CIs),
        offset(use_offset ? frame.offset : 0),
        use_diff_angle(use_diff_angle)
      {}

      void operator()() {
        try {
          const cart_t K = cart_t(0, 0, -parent.Kl);
          if (parent.mat_type == 3) { // 2-beam
            for (size_t i = 0; i < frame.strong_measured_beams.size(); i++) {
              size_t beam_idx = frame.strong_measured_beams[i];
              miller::index<> h = frame.indices[beam_idx];
              int ii = parent.mi_lookup.find_hkl(h);
              complex_t Fc = ii != -1 ? Fcs_k[ii] : 0;
              if (use_diff_angle) {
                frame.update_alpha(frame.beams[beam_idx].diffraction_angle);
              }
              complex_t ci = utils<FloatType>::calc_amp_2beam(
                h, parent.Fc2Ug * Fc,
                thickness,
                K, frame.RMf, frame.normal);
              Is[offset + i] = std::norm(ci);
              if (CIs.size() != 0) {
                CIs[offset + i] = ci;
              }
            }
            return;
          }
          // modified Ug
          if (parent.mat_type == 4) {
            af::versa<complex_t, af::mat_grid> A;
            af::shared<miller::index<> >
              s = af::select(frame.indices.const_ref(), frame.strong_beams.const_ref());

            utils<FloatType>::build_Ug_matrix(
              A, Fcs_k,
              parent.mi_lookup,
              s,
              parent.Fc2Ug
            );

            utils<FloatType>::modify_Ug_matrix(
              A, s, Fcs_k,
              parent.mi_lookup,
              af::select(frame.indices.const_ref(), frame.weak_beams.const_ref()),
              af::select(frame.excitation_errors.const_ref(), frame.weak_beams.const_ref()),
              parent.Fc2Ug);

            af::shared<FloatType> ExpDen;
            utils<FloatType>::build_eigen_matrix_modified(
              A, K, frame.normal,
              af::select(frame.gs.const_ref(), frame.strong_beams.const_ref()),
              af::select(frame.excitation_errors.const_ref(), frame.strong_beams.const_ref()),
              ExpDen);

            af::shared<complex_t> amps =
              utils<FloatType>::calc_amps_modified(A, ExpDen, thickness,
                frame.strong_measured_beams.size());
            
            for (size_t i = 0; i < amps.size(); i++) {
              Is[offset + i] = std::norm(amps[i]);
              if (CIs.size() != 0) {
                CIs[offset + i] = amps[i];
              }
            }
            return;
          }
          af::versa<complex_t, af::mat_grid> A;
          af::shared<miller::index<> >
            strong_indices = af::select(frame.indices.const_ref(),
              frame.strong_beams.const_ref());
          utils<FloatType>::build_Ug_matrix(
            A, Fcs_k,
            parent.mi_lookup,
            strong_indices,
            parent.Fc2Ug
          );
          af::shared<complex_t> amps;
          // 2013
          if (parent.mat_type == 0) {
            af::shared<FloatType> ExpDen;
            utils<FloatType>::build_eigen_matrix_2013(A, strong_indices, K,
              frame.RMf, frame.normal, ExpDen, parent.Fc2Ug
            );
            amps = utils<FloatType>::calc_amps_2013(A, ExpDen,
              thickness, frame.strong_measured_beams.size());
          }
          else if (parent.mat_type == 1) { // 2015
            af::shared<FloatType> M;
            utils<FloatType>::build_eigen_matrix_2015(A, strong_indices, K,
              frame.RMf, frame.normal, M, parent.Fc2Ug
            );
            const FloatType Kn = -frame.normal[2] * parent.Kl;
            amps = utils<FloatType>::calc_amps_2015(A, M,
              thickness, Kn, frame.strong_measured_beams.size());
          }
          else if (parent.mat_type == 2) { // ReciPro
            af::shared<FloatType> Pgs;
            utils<FloatType>::build_eigen_matrix_recipro(A, strong_indices, K,
              frame.RMf, frame.normal, Pgs, parent.Fc2Ug
            );
            amps = utils<FloatType>::calc_amps_recipro(A, Pgs, parent.Kvac,
              thickness, frame.strong_measured_beams.size());
          }
          for (size_t i = 0; i < amps.size(); i++) {
            Is[offset + i] = std::norm(amps[i]);
            if (CIs.size() != 0) {
              CIs[offset + i] = amps[i];
            }
          }
        }
        catch (smtbx::error const& e) {
          exception_.reset(new smtbx::error(e));
        }
        catch (std::exception const& e) {
          exception_.reset(new smtbx::error(e.what()));
        }
      }
      ed_n_shared_data const& parent;
      FrameInfo<FloatType>& frame;
      FloatType thickness;
      const af::shared<complex_t>& Fcs_k;
      af::shared<FloatType>& Is;
      af::shared<complex_t>& CIs;
      size_t offset;
      bool use_diff_angle;
      boost::scoped_ptr<smtbx::error> exception_;
    };

    af::shared<FloatType> process_frame_id(int frame_id,
      af::shared<complex_t> const& Fcs_k)
    {
      typename std::map<int, FrameInfo<FloatType>*>::const_iterator fi =
        frames_map.find(frame_id);
      SMTBX_ASSERT(fi != frames_map.end());
      FrameInfo<FloatType> &frame= *fi->second;
      af::shared<FloatType> Is_(frame.beams.size());
      af::shared<complex_t> CIs_;
      process_frame(*this, frame, Fcs_k, Is_, CIs_,
        false, false)(); // no offset - write to Is_ at 0;
      return Is_;
    }

    af::shared<PeakProfilePoint<FloatType> > build_profile(int frame_id,
      af::shared<complex_t> const& Fcs_k, FloatType angle,
      FloatType step)
    {
      typename std::map<int, FrameInfo<FloatType>*>::const_iterator fi =
        frames_map.find(frame_id);
      SMTBX_ASSERT(fi != frames_map.end());
      FrameInfo<FloatType>& frame = *fi->second;
      af::shared<FloatType> Is_(frame.beams.size());
      af::shared<complex_t> CIs_;
      af::shared<PeakProfilePoint<FloatType> > rv;
      FloatType alpha = frame.alpha;
      const cart_t K = cart_t(0, 0, -Kl);
      for (FloatType s = -angle; s <= angle; s+= step) {
        frame.update_alpha(alpha + s);
        process_frame(*this, frame, Fcs_k, Is_, CIs_,
          false, false)(); // no offset - write to Is_ at 0;
        for (size_t i = 0; i < frame.strong_measured_beams.size(); i++) {
          size_t beam_idx = frame.strong_measured_beams[i];
          miller::index<> h = frame.indices[beam_idx];
          cart_t g = frame.RMf * cart_t(h[0], h[1], h[2]);
          FloatType Sg = utils<FloatType>::calc_Sg(g, K);
          rv.push_back(PeakProfilePoint<FloatType>(Is_[i], Sg, frame.alpha));
        }
      }
      frame.update_alpha(alpha);
      return rv;
    }

    void process_frames_mt(af::shared<FloatType> &Is_,
      af::shared<complex_t>& CIs_,
      af::shared<complex_t> const& Fcs_k)
    {
      if (thread_n < 0) {
        thread_n = builder_base<FloatType>::get_available_threads();
      }
      boost::thread_group pool;
      typedef boost::shared_ptr<process_frame> frame_processor_t;
      size_t to = 0;
      for (size_t fi = 0; fi < frames.size(); fi += thread_n) {
        size_t t_end = std::min(thread_n, (int)(frames.size() - fi));
        if (t_end == 0) {
          break;
        }
        std::vector<frame_processor_t> accumulators;
        for (int thread_idx = 0; thread_idx < t_end; thread_idx++) {
          frame_processor_t pf(
            new process_frame(*this,
              frames[to],
              Fcs_k,
              Is_,
              CIs_, true, use_diff_angle)
          );
          accumulators.push_back(pf);
          pool.create_thread(boost::ref(*pf));
          to++;
        }
        pool.join_all();
        for (int thread_idx = 0; thread_idx < t_end; thread_idx++) {
          if (accumulators[thread_idx]->exception_) {
            throw* accumulators[thread_idx]->exception_.get();
          }
        }
      }
    }

    void build() {
      // Generate Fcs at current position
      {
        Is.resize(beam_n);
        Fcs.resize(beam_n);
        if (Fcs_k.size() != indices.size()) {
          Fcs_k.resize(indices.size());
          for (size_t ih = 0; ih < indices.size(); ih++) {
            Fcs_k[ih] = calc_one_h(indices[ih]);
          }
          // expand uniq Fc to frame indices
          size_t offset = 0;
          for (size_t i = 0; i < frames.size(); i++) {
            const af::shared<miller::index<> >& fidx = frames[i].indices;
            size_t measured = frames[i].strong_measured_beams.size();
            for (size_t i = 0; i < measured; i++) {
              long idx = mi_lookup.find_hkl(fidx[i]);
              Fcs[offset + i] = Fcs_k[idx];
            }
            offset += measured;
          }
        }
        // replacing dummy with Fc will allow to collect complex amplitudes
        af::shared<complex_t> dummy;
        process_frames_mt(Is, dummy, Fcs_k);
        if (!compute_grad) {
          return;
        }
      }
      af::shared<parameter*> params = reparamn.independent();
      af::shared<asu_parameter*> p_owners = reparamn.independent_owners(params);
        size_t param_n = 0;
      for (size_t i = 0; i < params.size(); i++) {
        param_n += params[i]->components().size();
      }
      design_matrix.resize(af::c_grid<2>(beam_n, param_n));
      // Generate Arrays for storing positive and negative epsilon Fcs for numerical
      // generation of gradients
      af::shared<complex_t> Fc_eps(indices.size()), dummy;
      af::shared<FloatType> Is_p(beam_n), Is_m(beam_n);
      FloatType t_eps = 2 * eps;
      for (size_t i = 0, n = 0; i < params.size(); i++) {
        parameter* param = params[i];
        af::ref<double> x = param->components();
        asu_parameter* cp = p_owners[i];
        for (size_t j = 0; j < x.size(); j++, n++) {
          x[j] += eps;
          if (cp != 0) {
            cp->evaluate(reparamn.unit_cell());
            cp->store(reparamn.unit_cell());
          }
          for (size_t i_h = 0; i_h < indices.size(); i_h++) {
            Fc_eps[i_h] = calc_one_h(indices[i_h]);
          }
          //Generate Fcs at x+eps
          process_frames_mt(Is_p, dummy, Fc_eps);

          x[j] -= t_eps;
          if (cp != 0) {
            cp->evaluate(reparamn.unit_cell());
            cp->store(reparamn.unit_cell());
          }
          for (size_t i_h = 0; i_h < indices.size(); i_h++) {
            Fc_eps[i_h] = calc_one_h(indices[i_h]);
          }
          //Generate Fcs at x-eps
          process_frames_mt(Is_m, dummy, Fc_eps);
          
          // compute grads
          for (size_t bi = 0; bi < beam_n; bi++) {
            FloatType grad_I = (Is_p[bi] - Is_m[bi]) / t_eps;
            design_matrix(bi, n) = grad_I;
          }

          // reset to original value
          x[j] += eps;
          if (cp != 0) {
            cp->store(reparamn.unit_cell());
          }
        }
      }
    }
  

    /* computes the position of the given miller index of the given
    frame in the uniform arrays
    */
    size_t find_hkl(int frame_id, miller::index<> const& h) const {
      typename std::map<int, lookup_ptr_t>::const_iterator i =
        frame_lookups.find(frame_id);
      SMTBX_ASSERT(i != frame_lookups.end());
      typename std::map<int, FrameInfo<FloatType> *>::const_iterator fi =
        frames_map.find(frame_id);
      long hi = i->second->find_hkl(h);
      SMTBX_ASSERT(hi >= 0);
      return hi + fi->second->offset;
    }

    reparametrisation const& reparamn;
    f_calc_function_base_t& f_calc_function;
    cctbx::xray::fc_correction<FloatType> const& fc_cr;
    sgtbx::space_group const& space_group;
    af::shared<miller::index<> > indices;
    FloatType Kvac, Kl, Fc2Ug, eps;
    int mat_type;
    size_t beam_n;
    /* to lookup an index in particular frame, have to keep a copy of the
    indices
    */
    typename std::map<int, lookup_ptr_t> frame_lookups;
    typename std::map<int, FrameInfo<FloatType>*> frames_map;
    af::shared<FrameInfo<FloatType> > frames;
    cctbx::xray::thickness<FloatType> const& thickness;
    bool compute_grad;
    // newly-calculated, aligned by frames
    af::shared<complex_t> Fcs, Fcs_k;
    af::shared<FloatType> Is;
    // 
    af::versa<FloatType, af::c_grid<2> > design_matrix;
    lookup_t mi_lookup;
    int thread_n;
    bool use_diff_angle;
  };

  template <typename FloatType>
  class f_calc_function_ed_n : public f_calc_function_base<FloatType> {
  public:
    typedef f_calc_function_base<FloatType> f_calc_function_base_t;
    typedef ed_n_shared_data<FloatType> data_t;
    typedef scitbx::vec3<FloatType> cart_t;

    f_calc_function_ed_n(ed_n_shared_data<FloatType> const& data)
      : data(data),
      index(~0)
    {}
    f_calc_function_ed_n(f_calc_function_ed_n const& other)
      : data(other.data),
      index(~0)
    {}

    virtual void compute(
      miller::index<> const& h,
      boost::optional<std::complex<FloatType> > const& f_mask = boost::none,
      twin_fraction<FloatType> const* fraction = 0,
      bool compute_grad = true)
    {
      SMTBX_ASSERT(fraction->tag >= 0);
      index = data.find_hkl(fraction->tag, h);
    }

    virtual boost::shared_ptr<f_calc_function_base_t> fork() const {
      return boost::shared_ptr<f_calc_function_base_t>(
        new f_calc_function_ed_n(*this));
    }

    virtual FloatType get_observable() const {
      return data.Is[index];
    }
    virtual std::complex<FloatType> get_f_calc() const {
      return data.Fcs[index];
    }
    virtual af::const_ref<FloatType> get_grad_observable() const {
      typedef af::versa_plain<FloatType> one_dim_type;
      typedef typename one_dim_type::accessor_type one_dim_accessor_type;
      one_dim_accessor_type a(data.design_matrix.accessor().n_columns());
      return af::const_ref<FloatType>(&data.design_matrix(index, 0), a);
    }

    virtual bool raw_gradients() const { return false; }
  private:
    ed_n_shared_data<FloatType> const& data;
    size_t index;
  };

}}}

#endif // GUARD
