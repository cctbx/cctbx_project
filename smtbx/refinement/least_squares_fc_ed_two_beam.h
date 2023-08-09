#ifndef SMTBX_REFINEMENT_LEAST_SQUARES_FC_ED_TB_H
#define SMTBX_REFINEMENT_LEAST_SQUARES_FC_ED_TB_H

#include <cctbx/miller/lookup_utils.h>
#include <cctbx/xray/thickness.h>

#include <smtbx/refinement/least_squares_fc.h>
#include <smtbx/import_scitbx_af.h>
#include <scitbx/vec3.h>
#include <scitbx/constants.h>
#include <smtbx/ED/ed_data.h>
#include <smtbx/ED/kinematic.h>

namespace smtbx { namespace refinement {
namespace least_squares {
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
      // Kl, Fc2Ug, epsilon
      af::shared<FloatType> const& params,
      bool compute_grad,
      bool do_build = true)
      : Jt_matching_grad_fc(Jt_matching_grad_fc),
      f_calc_function(f_calc_function),
      space_group(space_group),
      Kl(params[1]),
      Fc2Ug(params[2]),
      frames(frames),
      thickness(thickness),
      compute_grad(compute_grad),
      thread_n(static_cast<int>(params[5]))
    {
      K = cart_t(0, 0, -Kl);
      // build lookups for each frame + collect all indices and they diffs
      af::shared<miller::index<> > all_indices;
      size_t offset = 0;
      // treat equivalents independently inside the frames
      sgtbx::space_group P1("P 1");
      for (size_t i = 0; i < frames.size(); i++) {
        FrameInfo<FloatType>& frame = frames[i];
        frames_map.insert(std::make_pair(frame.id, &frame));
        af::shared<miller::index<> > indices =
          af::select(frame.indices.const_ref(),
            frame.strong_measured_beams.const_ref());
        lookup_ptr_t mi_l(new lookup_t( indices.const_ref(), P1, true));
        frame_lookups.insert(std::make_pair(frame.id, mi_l));
        all_indices.extend(indices.begin(), indices.end());
      }
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

    const scitbx::sparse::matrix<FloatType>& Jt_matching_grad_fc;
    f_calc_function_base_t& f_calc_function;
    sgtbx::space_group const& space_group;
    af::shared<miller::index<> > indices;
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

    f_calc_function_ed_two_beam(f_calc_function_ed_two_beam const & other)
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
      Fsq = 0;
      index = data.mi_lookup.find_hkl(h);
      if (index == -1) {
        if (!data.space_group.is_sys_absent(h)) {
          SMTBX_ASSERT(index >= 0)(h.as_string());
        }
        Fc = 0;
        observable_updated = true;
      }
      else {
        observable_updated = false;
        Fc = data.Fcs_kin[index];
      }
      typename std::map<int, FrameInfo<FloatType>*>::const_iterator fi =
        data.frames_map.find(fraction->tag);
      SMTBX_ASSERT(fi != data.frames_map.end());
      frame = fi->second;
      this->h = h;
      if (compute_grad) {
        grads.resize(data.design_matrix_kin.accessor().n_columns() +
        (data.thickness.grad ? 1 : 0));
        std::fill(grads.begin(), grads.end(), 0);
      }
      else {
        grads.resize(0);
      }
      computed = true;
    }

    virtual boost::shared_ptr<f_calc_function_base_t> fork() const {
      return boost::shared_ptr<f_calc_function_base_t>(
        new f_calc_function_ed_two_beam(*this));
    }

    FloatType get_observable_int() const {
      FloatType angle = scitbx::deg_as_rad(1.0),
        step = scitbx::deg_as_rad(0.075);
      int steps = round(angle / step);
      FloatType I1 = -1, g1 = -1, I_sum = 0;
      af::shared<complex_t> grads1, grads2;
      FloatType da = frame->get_diffraction_angle(h, data.Kl);
      for (int st = -steps; st <= steps; st++) {
        std::pair<mat3_t, cart_t> r = frame->compute_RMf_N(
          da + st * step);
        cart_t K_g = r.first * cart_t(h[0], h[1], h[2]) + data.K;
        FloatType I;
        if (!grads.empty()) {
          size_t n_param = data.design_matrix_kin.accessor().n_columns();
          if (n_param > 0) {
            grads1.assign(&data.design_matrix_kin(index, 0),
              &data.design_matrix_kin(index, n_param));
          }
          else {
            grads1.clear();
          }
          I = std::norm(utils<FloatType>::calc_2beam_ext(
            h, Fc, grads1,
            data.thickness.value,
            data.K, r.first, r.second, data.thickness.grad));
        }
        else {
          I = std::norm(utils<FloatType>::calc_amp_2beam(
            h, Fc,
            data.thickness.value,
            data.K, r.first, r.second));
        }
        FloatType g = K_g.length();
        if (g1 >= 0) {
          FloatType st = std::abs(g - g1) / 2;
          I_sum += (I + I1) * st;
          for (size_t i = 0; i < grads.size(); i++) {
            grads[i] += (grads1[i].real() + grads2[i].real()) * st;
          }
        }
        g1 = g;
        I1 = I;
        if (!grads.empty()) {
          grads2 = grads1.deep_copy();
        }
      }
      return I_sum;
    }

    virtual FloatType get_observable() const {
      SMTBX_ASSERT(frame != 0);
      Fsq = get_observable_int();
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
    mutable bool observable_updated, computed;
      mutable complex_t Fc;
    mutable FloatType Fsq;
    mutable af::shared<FloatType> grads;
    miller::index<> h;
  };

}}}

#endif // GUARD
