#ifndef SMTBX_REFINEMENT_LEAST_SQUARES_FC_ED_TB_H
#define SMTBX_REFINEMENT_LEAST_SQUARES_FC_ED_TB_H

#include <cctbx/miller/lookup_utils.h>
#include <cctbx/xray/thickness.h>

#include <smtbx/refinement/least_squares_fc.h>
#include <scitbx/vec3.h>
#include <scitbx/constants.h>
#include <smtbx/ED/n_beam.h>

namespace smtbx {  namespace refinement  { namespace least_squares
{
  using namespace smtbx::ED;

  template <typename FloatType>
  struct two_beam_shared_data {
    ED_UTIL_TYPEDEFS;
    typedef f_calc_function_base<FloatType> f_calc_function_base_t;
    typedef boost::shared_ptr< f_calc_function_base_t> f_calc_function_base_ptr_t;
    typedef builder_base<FloatType> data_t;
    typedef af::shared<const BeamInfo<FloatType>*> beam_at;
    typedef std::pair<int, af::shared<const BeamInfo<FloatType>*> > beam_me;
    typedef boost::shared_ptr<lookup_t> lookup_ptr_t;
    typedef cctbx::xray::fc_correction<FloatType> fc_correction_t;
    typedef boost::shared_ptr< fc_correction_t> fc_correction_ptr_t;

    two_beam_shared_data(const scitbx::sparse::matrix<FloatType>&
      Jt_matching_grad_fc,
      f_calc_function_base_t& f_calc_function,
      sgtbx::space_group const& space_group,
      bool anomalous_flag,
      af::shared<FrameInfo<FloatType> > frames,
      cctbx::xray::thickness<FloatType> const& thickness,
      const RefinementParams<FloatType>& params,
      bool compute_grad,
      bool do_build = true)
      : Jt_matching_grad_fc(Jt_matching_grad_fc),
      f_calc_function(f_calc_function),
      space_group(space_group),
      params(params),
      Kl(params.getKl()),
      Fc2Ug(params.getFc2Ug()),
      frames(frames),
      thickness(thickness),
      compute_grad(compute_grad),
      thread_n(params.getThreadN())
    {
      K = cart_t(0, 0, -Kl);
      // build lookups for each frame + collect all indices and they diffs
      af::shared<miller::index<> > all_indices;
      // treat equivalents independently inside the frames
      sgtbx::space_group P1("P 1");
      if (params.getBeamN() == 2) {
        for (size_t i = 0; i < frames.size(); i++) {
          FrameInfo<FloatType>& frame = frames[i];
          frames_map.insert(std::make_pair(frame.id, &frame));
          af::shared<miller::index<> > indices =
            af::select(frame.indices.const_ref(),
              frame.strong_measured_beams.const_ref());
          lookup_ptr_t mi_l(new lookup_t(indices.const_ref(), P1, true));
          frame_lookups.insert(std::make_pair(frame.id, mi_l));
          all_indices.extend(indices.begin(), indices.end());
        }
      }
      else {
        for (size_t i = 0; i < frames.size(); i++) {
          FrameInfo<FloatType>& frame = frames[i];
          frames_map.insert(std::make_pair(frame.id, &frame));
          for (size_t hi = 0; hi < frame.strong_beams.size(); hi++) {
            const miller::index<>& h = frame.indices[frame.strong_beams[hi]];
            all_indices.push_back(h);
            all_indices.push_back(-h);
            for (size_t hj = hi + 1; hj < frame.strong_beams.size(); hj++) {
              const miller::index<>& k = frame.indices[frame.strong_beams[hj]];
              all_indices.push_back(h - k);
              all_indices.push_back(k - h);
            }
          }
          lookup_ptr_t mi_l(new lookup_t(
            af::select(frame.indices.const_ref(),
              frame.strong_measured_beams.const_ref()).const_ref(),
            P1,
            true));
          frame_lookups.insert(std::make_pair(frames[i].id, mi_l));
        }
      }
      // a tricky way of getting unique only...
      mi_lookup = lookup_t(
        all_indices.const_ref(),
        P1,
        anomalous_flag);
      indices = mi_lookup.get_unique();
      mi_lookup = lookup_t(
        indices.const_ref(),
        P1,
        anomalous_flag);
      if (do_build) {
        build();
      }
    }

    ~two_beam_shared_data() {
    }

    void do_build_kin_mt() {
      if (thread_n < 0) {
        thread_n = builder_base<FloatType>::get_available_threads();
      }
      build_kin_mt(thread_n, Jt_matching_grad_fc, Fc2Ug, f_calc_function,
        indices, Fcs_kin, design_matrix_kin, compute_grad);
    }

    void build() {
      if (Fcs_kin.size() != indices.size()) {
        if (compute_grad) {
          size_t cn = Jt_matching_grad_fc.n_rows() - (thickness.grad ? 1 : 0);
          if (cn > 0) {
            design_matrix_kin.resize(
              af::mat_grid(indices.size(), cn));
          }
        }
        Fcs_kin.resize(indices.size());
        do_build_kin_mt();
      }
    }

    scitbx::sparse::matrix<FloatType> Jt_matching_grad_fc;
    f_calc_function_base_t& f_calc_function;
    sgtbx::space_group const& space_group;
    af::shared<miller::index<> > indices;
    RefinementParams<FloatType> params;
    FloatType Kl, Fc2Ug;
    cart_t K;
    /* to lookup an index in particular frame, have to keep a copy of the
    indices
    */
    typename std::map<int, lookup_ptr_t> frame_lookups;
    typename std::map<int, FrameInfo<FloatType>*> frames_map;
    af::shared<FrameInfo<FloatType> > frames;
    cctbx::xray::thickness<FloatType> const& thickness;
    bool compute_grad;
    af::shared<complex_t> Fcs_kin;
    // 
    cmat_t design_matrix_kin;
    lookup_t mi_lookup;
    int thread_n;
  };

  template <typename FloatType>
  class f_calc_function_ed_two_beam : public f_calc_function_base<FloatType> {
  public:
    ED_UTIL_TYPEDEFS;
    typedef f_calc_function_base<FloatType> f_calc_function_base_t;
    typedef two_beam_shared_data<FloatType> data_t;
    typedef af::versa_plain<FloatType> one_dim_type;
    typedef typename one_dim_type::accessor_type one_dim_accessor_type;

    f_calc_function_ed_two_beam(data_t const& data)
      : data(data),
      index(-1),
      observable_updated(false),
      computed(false)
    {
      frame = 0;
    }

    f_calc_function_ed_two_beam(f_calc_function_ed_two_beam const& other)
      : data(other.data),
      observable_updated(false),
      computed(false)
    {}

    virtual void compute(
      miller::index<> const& h,
      boost::optional<complex_t> const& f_mask = boost::none,
      twin_fraction<FloatType> const* fraction = 0,
      bool compute_grad = true)
    {
      SMTBX_ASSERT(fraction != 0);
      index = data.mi_lookup.find_hkl(h);
      if (index == -1) {
        if (!data.space_group.is_sys_absent(h)) {
          SMTBX_ASSERT(index >= 0)(h.as_string());
        }
        Fc = 0;
        Fsq = 0;
        observable_updated = true;
      }
      else {
        observable_updated = false;
        Fc = data.Fcs_kin[index];
        Fsq = std::norm(Fc);
      }
      typename std::map<int, FrameInfo<FloatType>*>::const_iterator fi =
        data.frames_map.find(fraction->tag);
      SMTBX_ASSERT(fi != data.frames_map.end());
      frame = fi->second;
      this->h = h;
      this->compute_grad = compute_grad;
      computed = true;
    }

    virtual boost::shared_ptr<f_calc_function_base_t> fork() const {
      return boost::shared_ptr<f_calc_function_base_t>(
        new f_calc_function_ed_two_beam(*this));
    }

    // for derivatives testing
    FloatType calc(FloatType Fsq, cart_t const& K, FloatType t,
      FloatType Sg, cart_t const& n) const
    {
      FloatType Kl = K.length(),
        t_part = ((scitbx::constants::pi * t) / (K * n));
      FloatType X = scitbx::fn::pow2(Kl * Sg) + Fsq,
        sin_part = scitbx::fn::pow2(std::sin(t_part * std::sqrt(X)));
      return Fsq * sin_part / X;
    }

    // Acta Cryst. (2013). A69, 171–188
    FloatType get_observable_pltns_2013() const {
      FloatType da = frame->get_diffraction_angle(h, data.Kl);
      size_t n_param = data.Jt_matching_grad_fc.n_rows();
      size_t coln = data.design_matrix_kin.accessor().n_columns();
      af::shared<FloatType> grads_sum(coln);
      FloatType I1 = -1, g1 = -1, grad_fsq1 = 0, grad_t1 = 0, I_sum = 0,
        dT_sum = 0;
    af::shared<FloatType> angles = frame->get_angles(da,
      data.params.getIntSpan(),
      data.params.getIntStep());
    for (size_t ai = 0; ai < angles.size(); ai++) {
        std::pair<mat3_t, cart_t> r = frame->compute_RMf_N(angles[ai]);
        cart_t K_g = r.first * cart_t(h[0], h[1], h[2]) + data.K;
        FloatType g = K_g.length();
        af::shared<FloatType> res = get_observable_pltns_2013_(angles[ai]);
        if (g1 >= 0) {
          FloatType d = std::abs(g - g1) / 2;
          I_sum += (res[0] + I1) * d;
          if (compute_grad) {
            for (size_t i = 0; i < coln; i++) {
              grads_sum[i] += (res[1] + grad_fsq1) * d;
            }
            if (data.thickness.grad) {
              dT_sum += (res[2] + grad_t1) * d;
            }
          }
        }
        I1 = res[0];
        g1 = g;
        grad_fsq1 = res[1];
        grad_t1 = res[2];
      }
      if (compute_grad) {
        grads = af::shared<FloatType>(n_param);
        for (size_t i = 0; i < coln; i++) {
          complex_t dFc = data.design_matrix_kin(index, i);
          FloatType dFc_sq = 2.0 * (dFc.real() * Fc.real() + dFc.imag() * Fc.imag());
          grads[i] = dFc_sq * grads_sum[i];
        }
        if (data.thickness.grad) {
          int grad_index = data.thickness.grad_index;
          SMTBX_ASSERT(!(grad_index < 0 || grad_index >= grads.size()));
          grads[grad_index] = dT_sum;
        }

      }
      return I_sum;
    }
    // returns I, d_I_d_F_sq, d_I_d_EDT
    af::shared<FloatType> get_observable_pltns_2013_(FloatType angle) const {
      FloatType Kl = data.Kl;
      std::pair<mat3_t, cart_t> FI = frame->compute_RMf_N(angle);

      cart_t g = FI.first * cart_t(h[0], h[1], h[2]);
      cart_t K = cart_t(0, 0, -Kl);

      FloatType Sg = (Kl * Kl - (K + g).length_sq()) / (2 * Kl);
      FloatType X = scitbx::fn::pow2(Kl * Sg) + Fsq,
        t = data.thickness.value,
        t_part = ((scitbx::constants::pi * t) / (K * FI.second)),
        P = t_part * std::sqrt(X);
      FloatType sin_part = scitbx::fn::pow2(std::sin(P)),
        I = Fsq * sin_part / X;
      af::shared<FloatType> rv;
      rv.push_back(I);
      // update gradients...
      if (compute_grad && index >= 0) {
        if (data.Jt_matching_grad_fc.n_cols() > 0) {
          FloatType grad_fsq = sin_part / X + Fsq *
            (2 * std::sin(P) * std::cos(P) * (t_part * std::pow(X, -0.5) / 2.0) * X - sin_part) / (X * X);
          rv.push_back(grad_fsq);
        }
        else {
          rv.push_back(0);
        }
        if (data.thickness.grad) {
          rv.push_back(Fsq * 2 * std::sin(P) * std::cos(P) * P / t / X);
        }
        else {
          rv.push_back(0);
        }
        /* Testing derivatives shows consistensy with analytical expressions above */
        //FloatType eps = 1e-6;
        //FloatType v1 = calc(Fsq - eps, K, t, Sg, frame->normal);
        //FloatType v2 = calc(Fsq + eps, K, t, Sg, frame->normal);
        //FloatType diff = (v2 - v1) / (2*eps);
        //
        //v1 = calc(Fsq, K, t - eps, Sg, frame->normal);
        //v2 = calc(Fsq, K, t + eps, Sg, frame->normal);
        //diff = (v2 - v1) / (2 * eps);
        //FloatType grad_t = Fsq * 2 * std::sin(P) * std::cos(P) * P / thickness.value / X;
        //FloatType v = calc(Fsq, K, Sg, t, frame->normal);
        //v = 0; // allow for a breakpoint here
      }
      return rv;
    }

    FloatType get_observable_N() const {
      FloatType da = frame->get_diffraction_angle(h, data.Kl);
      af::shared<FloatType> angles;
      if (data.params.isAngleInt()) {
        angles = frame->get_angles(da,
          data.params.getIntSpan(),
          data.params.getIntStep());
      }
      else {
        angles = frame->get_angles_Sg(h, data.Kl,
          data.params.getIntSpan(),
          data.params.getIntStep());
      }

      dyn_calculator_n_beam<FloatType> n_beam_dc(data.params.getBeamN(),
        data.params.getMatrixType(),
        *frame, data.K, data.thickness.value,
        data.params.useNBeamSg(), data.params.getNBeamWght());

      n_beam_dc.init(h, da, data.Fcs_kin, data.mi_lookup);

      af::shared<cmat_t> Ds_kin;
      mat_t D_dyn;
      if (compute_grad) {
        utils<FloatType>::build_D_matrices(data.mi_lookup, n_beam_dc.indices,
          data.design_matrix_kin, Ds_kin);
      }
      size_t n_param = data.Jt_matching_grad_fc.n_rows();
      size_t coln = data.design_matrix_kin.accessor().n_columns() +
        (data.thickness.grad ? 1 : 0);
      af::shared<FloatType> grads_sum(coln),
        grads1;
      FloatType I1 = -1, g1 = -1, I_sum = 0;
      for (size_t ai=0; ai < angles.size(); ai++) {
        std::pair<mat3_t, cart_t> r = frame->compute_RMf_N(angles[ai]);
        cart_t g = r.first * cart_t(h[0], h[1], h[2]);
        cart_t K_g = g + data.K;
        FloatType K_g_l = K_g.length();

        FloatType I;
        if (compute_grad) {
          I = std::norm(
            n_beam_dc.calc_amp_ext(r, Ds_kin, data.thickness.grad, D_dyn));
        }
        else {
          I = std::norm(n_beam_dc.calc_amp(r));
        }
        if (g1 >= 0) {
          FloatType d = std::abs(K_g_l - g1) / 2;
          if (d == 0) {
            continue;
          }
          I_sum += (I + I1) * d;
          if (compute_grad) {
            for (size_t i = 0; i < coln; i++) {
              grads_sum[i] += (grads1[i] + D_dyn(0, i)) * d;
            }
          }
        }
        I1 = I;
        g1 = K_g_l;
        if (compute_grad) {
          grads1 = af::shared<FloatType>(&D_dyn(0, 0), &D_dyn(0, n_param));
        }
      }
      if (compute_grad) {
        grads = grads_sum;
      }
      return I_sum;
    }

    virtual FloatType get_observable() const {
      if (observable_updated || !computed) {
        return Fsq;
      }
      SMTBX_ASSERT(frame != 0);
      if (data.params.getBeamN() == 2) {
        Fsq = get_observable_pltns_2013();
      }
      else {
        Fsq = get_observable_N();
      }
      observable_updated = true;
      return Fsq;
    }

    virtual complex_t get_f_calc() const {
      if (!observable_updated) {
        get_observable();
      }
      return Fc;
    }
    virtual af::const_ref<complex_t> get_grad_f_calc() const {
      SMTBX_NOT_IMPLEMENTED();
      throw 1;
    }

    virtual af::const_ref<FloatType> get_grad_observable() const {
      if (!observable_updated) {
        get_observable();
      }
      return grads.const_ref();
    }

    virtual bool raw_gradients() const { return false; }
  private:
    data_t const& data;
    long index;
    const FrameInfo<FloatType>* frame;
    bool compute_grad;
    mutable bool observable_updated, computed;
    mutable complex_t Fc;
    mutable FloatType Fsq;
    mutable af::shared<FloatType> grads;
    miller::index<> h;
  };

}}}

#endif // GUARD
