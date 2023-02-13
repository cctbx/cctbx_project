#ifndef SMTBX_REFINEMENT_LEAST_SQUARES_FC_ED_TB_H
#define SMTBX_REFINEMENT_LEAST_SQUARES_FC_ED_TB_H

#include <cctbx/miller/lookup_utils.h>
#include <cctbx/xray/thickness.h>

#include <smtbx/refinement/least_squares_fc.h>
#include <smtbx/import_scitbx_af.h>
#include <scitbx/vec3.h>
#include <scitbx/constants.h>
#include <smtbx/ED/ed_data.h>

namespace smtbx { namespace refinement {
namespace least_squares {
  using namespace smtbx::ED;

  template <typename FloatType>
  class f_calc_function_ed_two_beam : public f_calc_function_base<FloatType> {
  public:
    typedef f_calc_function_base<FloatType> f_calc_function_base_t;
    typedef scitbx::vec3<FloatType> cart_t;
    typedef builder_base<FloatType> data_t;

    f_calc_function_ed_two_beam(data_t const& data,
      sgtbx::space_group const& space_group,
      bool anomalous_flag,
      af::shared<FrameInfo<FloatType> > const& frames,
      cctbx::xray::thickness<FloatType> const& thickness,
      af::shared<FloatType> const& params)
      : data(data),
      space_group(space_group),
      frames(frames),
      Kl(params[0]),
      Fc2Ug(params[1]),
      thickness(thickness),
      index(-1),
      observable_updated(false),
      computed(false)
    {
      f_calc = data.f_calc();
      observables = data.observables();
      design_matrix = data.design_matrix();
      mi_lookup = miller::lookup_utils::lookup_tensor<FloatType>(
        data.reflections().indices().const_ref(),
        space_group,
        anomalous_flag);
      for (size_t i = 0; i < frames.size(); i++) {
        frames_map.insert(std::make_pair(frames[i].id, &frames[i]));
      }
      frame = 0;
    } 
    f_calc_function_ed_two_beam(f_calc_function_ed_two_beam const & other)
      : data(other.data),
      space_group(other.space_group),
      Kl(other.Kl),
      Fc2Ug(other.Fc2Ug),
      frames(other.frames),
      frames_map(other.frames_map),
      thickness(other.thickness),
      mi_lookup(other.mi_lookup),
      observable_updated(false),
      computed(false)
    {
      f_calc = data.f_calc();
      observables = data.observables();
      design_matrix = data.design_matrix();
    }

    virtual void compute(
      miller::index<> const& h,
      boost::optional<std::complex<FloatType> > const& f_mask = boost::none,
      twin_fraction<FloatType> const* fraction = 0,
      bool compute_grad = true)
    {
      SMTBX_ASSERT(fraction != 0);
      index = mi_lookup.find_hkl(h);
      if (index == -1) {
        if (!space_group.is_sys_absent(h)) {
          SMTBX_ASSERT(index >= 0)(h.as_string());
        }
        Fc = 0;
        Fsq = 0;
        observable_updated = true;
      }
      else {
        observable_updated = false;
        Fc = f_calc[index];
        Fsq = observables[index] * Fc2Ug * Fc2Ug;
      }
      typename std::map<int, const FrameInfo<FloatType>*>::const_iterator fi =
        frames_map.find(fraction->tag);
      SMTBX_ASSERT(fi != frames_map.end());
      frame = fi->second;
      this->h = h;
      ratio = 1;
      if (compute_grad) {
        grads.resize(design_matrix.accessor().n_columns());
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

    /* IT Vol C, Chapter 4.3 */
    FloatType calc_test(FloatType U, FloatType t, FloatType Sg) const {
      FloatType sig = 0.0007289;
      FloatType delta_k = U * sig / scitbx::constants::pi;
      FloatType s_eff = std::sqrt(Sg * Sg + delta_k * delta_k);
      FloatType x1 = scitbx::fn::pow2(sig * U * sin(scitbx::constants::pi * t * s_eff));
      x1 /= scitbx::fn::pow2(scitbx::constants::pi * s_eff);
      return x1;
    }

    virtual FloatType get_observable() const {
      if (observable_updated || !computed) {
        return Fsq;
      }
      SMTBX_ASSERT(frame != 0);
      int coln = design_matrix.accessor().n_columns();
      cart_t g = frame->RMf * cart_t(h[0], h[1], h[2]);
      cart_t K = cart_t(0, 0, -Kl);
      FloatType Sg = (Kl * Kl - (K + g).length_sq()) / (2 * Kl);
      FloatType X = scitbx::fn::pow2(Kl * Sg) + Fsq,
        t = thickness.value,
        t_part = ((scitbx::constants::pi * t) / (K * frame->normal)),
        P = t_part * std::sqrt(X);
      FloatType sin_part = scitbx::fn::pow2(std::sin(P)),
        I = Fsq * sin_part / X;
      // update gradients...
      if (grads.size() > 0) {
        if (index >= 0) {
          FloatType grad_fsq = sin_part/X + Fsq *
            (2 * std::sin(P) * std::cos(P) * (t_part * std::pow(X, -0.5) / 2.0) * X - sin_part)/(X*X);
          for (int i = 0; i < coln; i++) {
            grads[i] = design_matrix(index, i) * grad_fsq;
          }
          if (thickness.grad) {
            int grad_index = thickness.grad_index;
            SMTBX_ASSERT(!(grad_index < 0 || grad_index >= grads.size()));
            grads[grad_index] = Fsq * 2 * std::sin(P) * std::cos(P) * P / thickness.value / X;
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
        else {
          std::fill(grads.begin(), grads.end(), 0);
        }
      }
      observable_updated = true;
      ratio = sin_part / X;
      return (Fsq = I);
    }
    virtual std::complex<FloatType> get_f_calc() const {
      if (!observable_updated) {
        get_observable();
      }
      return Fc;
    }
    // updated observable over original observable
    FloatType get_ratio() const {
      if (!observable_updated) {
        get_observable();
      }
      return ratio;
    }
    virtual af::const_ref<FloatType> get_grad_observable() const {
      if (!computed) {
        typedef af::versa_plain<FloatType> one_dim_type;
        typedef typename one_dim_type::accessor_type one_dim_accessor_type;
        one_dim_accessor_type a(design_matrix.accessor().n_columns());
        return af::const_ref<FloatType>(&design_matrix(index, 0), a);
        }
      if (!observable_updated) {
        get_observable();
      }
      return grads.const_ref();
    }

    virtual bool raw_gradients() const { return false; }
  private:
    data_t const& data;
    sgtbx::space_group const& space_group;
    af::shared<FrameInfo<FloatType> > frames;
    FloatType Kl, Fc2Ug;
    cctbx::xray::thickness<FloatType> const& thickness;
    af::shared<std::complex<FloatType> > f_calc;
    af::shared<FloatType> observables;
    af::shared<FloatType> weights;
    af::versa<FloatType, af::c_grid<2> > design_matrix;
    miller::lookup_utils::lookup_tensor<FloatType> mi_lookup;
    std::map<int, const FrameInfo<FloatType>*> frames_map;
    long index;
    const FrameInfo<FloatType>* frame;
    mutable bool observable_updated, computed;
    mutable std::complex<FloatType> Fc;
    mutable FloatType Fsq, ratio;
    mutable af::shared<FloatType> grads;
    miller::index<> h;
  };

}}}

#endif // GUARD
