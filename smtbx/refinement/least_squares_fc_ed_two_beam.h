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
      FloatType wavelength,
      bool anomalous_flag,
      scitbx::mat3<FloatType> const& UB,
      af::shared<FrameInfo<FloatType> > const& frames,
      af::shared<BeamInfo<FloatType> > const& beams,
      cctbx::xray::thickness<FloatType> const& thickness,
      double maxSg)
      : data(data),
      space_group(space_group),
      wavelength(wavelength),
      UB(UB),
      frames(frames),
      beams(beams),
      thickness(thickness),
      maxSg(maxSg),
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
    } 
    f_calc_function_ed_two_beam(f_calc_function_ed_two_beam const & other)
      : data(other.data),
      space_group(other.space_group),
      wavelength(other.wavelength),
      UB(other.UB),
      frames(other.frames),
      beams(other.beams),
      thickness(other.thickness),
      maxSg(other.maxSg),
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
      SMTBX_ASSERT(fraction != 0 &&
        fraction->tag >= 0 && fraction->tag < frames.size());
      index = mi_lookup.find_hkl(h);
      if (index == -1) {
        if (!space_group.is_sys_absent(h)) {
          SMTBX_ASSERT(index >= 0)(index);
        }
      }
      observable_updated = false;
      frame_index = fraction->tag;
      this->h = h;
      Fc = index < 0 ? 0 : f_calc[index];
      Fsq = index < 0 ? 0 : observables[index];
      ratio = 1;
      if (compute_grad) {
        grads.resize(design_matrix.accessor().n_columns()+
          (thickness.grad ? 1 : 0));
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
    FloatType calc(FloatType U, cart_t const& K, FloatType t,
      FloatType Sg, cart_t const& n) const
    {
      FloatType Kl = K.length(),
        t_part = ((scitbx::constants::pi * t) / (K * n)),
        Fsq = U*U;
      FloatType X = scitbx::fn::pow2(Kl * Sg) + Fsq,
        sin_part = scitbx::fn::pow2(std::sin(t_part * std::sqrt(X)));
      return Fsq * sin_part / X;
    }

    virtual FloatType get_observable() const {
      if (observable_updated || !computed) {
        return Fsq;
      }
      SMTBX_ASSERT(frame_index >= 0)(frame_index);
      int coln = design_matrix.accessor().n_columns();
      const FrameInfo<FloatType>& frame = frames[frame_index];
      cart_t g = frame.RM * UB * cart_t(h[0], h[1], h[2]);
      FloatType Kl = 1.0 / wavelength;
      cart_t K = cart_t(0, 0, -Kl);
      FloatType Sg = (Kl * Kl - scitbx::fn::pow2((K + g).length())) / (2 * Kl);
      if (std::abs(Sg) > maxSg) {
        observable_updated = true;
        if (grads.size() > 0) {
          if (index >= 0) {
            memcpy(grads.begin(), &design_matrix(index, 0), sizeof(FloatType)*coln);
            if (thickness.grad) {
              grads[thickness.grad_index] = 0;
            }
          }
          else {
            std::fill(grads.begin(), grads.end(), 0);
          }
        }
        return Fsq;
      }
      FloatType X = scitbx::fn::pow2(Kl * Sg) + Fsq,
        t = thickness.value,
        t_part = ((scitbx::constants::pi * t) / (K * frame.normal)),
        P = t_part * std::sqrt(X);
      FloatType sin_part = scitbx::fn::pow2(std::sin(P)),
        I = Fsq * sin_part / X;
      // update gradients...
      if (grads.size() > 0) {
        if (index >= 0) {
          FloatType U = std::sqrt(Fsq);
          FloatType grad_fc = 2*U*sin_part/X + Fsq *
            (2 * std::sin(P) * std::cos(P) * (t_part * std::pow(X, -0.5) * U) * X - sin_part * 2 * U)/(X*X);
          for (int i = 0; i < coln; i++) {
            grads[i] = design_matrix(index, i) * grad_fc;
          }
          if (thickness.grad) {
            int grad_index = thickness.grad_index;
            SMTBX_ASSERT(!(grad_index < 0 || grad_index >= grads.size()));
            grads[grad_index] = Fsq * 2 * std::sin(P) * std::cos(P) * P / thickness.value / X;
          }

          /* Testing derivatives shows consistensy with analytical expressions above
          FloatType eps = 1e-6;
          FloatType v1 = calc(U - eps, K, t, Sg, frame.normal);
          FloatType v2 = calc(U + eps, K, t, Sg, frame.normal);
          FloatType diff = (v2 - v1) / (2*eps);
          
          v1 = calc(U, K, t - eps, Sg, frame.normal);
          v2 = calc(U, K, t + eps, Sg, frame.normal);
          diff = (v2 - v1) / (2 * eps);
          FloatType grad_t = Fsq * 2 * std::sin(P) * std::cos(P) * P / thickness.value / X;
          FloatType v = calc(U, K, Sg, t, frame.normal);
          v = 0; // allow for a breakpoint here
          */
        }
        else {
          std::fill(grads.begin(), grads.end(), 0);
        }
      }
      observable_updated = true;
      ratio = sin_part / X;
      Fc *= std::sqrt(std::abs(ratio));
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
    FloatType wavelength;
    scitbx::mat3<FloatType> UB;
    af::shared<FrameInfo<FloatType> > frames;
    af::shared<BeamInfo<FloatType> > beams;
    cctbx::xray::thickness<FloatType> const& thickness;
    FloatType maxSg;
    af::shared<std::complex<FloatType> > f_calc;
    af::shared<FloatType> observables;
    af::shared<FloatType> weights;
    af::versa<FloatType, af::c_grid<2> > design_matrix;
    miller::lookup_utils::lookup_tensor<FloatType> mi_lookup;
    long index;
    int frame_index;
    mutable bool observable_updated, computed;
    mutable std::complex<FloatType> Fc;
    mutable FloatType Fsq, ratio;
    mutable af::shared<FloatType> grads;
    miller::index<> h;
  };

}}}

#endif // GUARD
