#ifndef SMTBX_REFINEMENT_LEAST_SQUARES_FC_ED_H
#define SMTBX_REFINEMENT_LEAST_SQUARES_FC_ED_H

#include <cctbx/xray/thickness.h>

#include <smtbx/refinement/least_squares_fc.h>
#include <smtbx/import_scitbx_af.h>
#include <scitbx/vec3.h>
#include <smtbx/ED/ed_data.h>

namespace smtbx {  namespace refinement {
namespace least_squares {
  using namespace smtbx::ED;

  template <typename FloatType>
  class f_calc_function_ed : public f_calc_function_base<FloatType> {
  public:
    typedef f_calc_function_base<FloatType> f_calc_function_base_t;
    typedef scitbx::vec3<FloatType> cart_t;
    typedef builder_base<FloatType> data_t;
          
    f_calc_function_ed(data_t const& data,
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
    f_calc_function_ed(f_calc_function_ed const& other)
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
        grads.resize(design_matrix.accessor().n_columns() +
          (thickness.grad ? 1 : 0));
      }
      else {
        grads.resize(0);
      }
      computed = true;
    }

    virtual boost::shared_ptr<f_calc_function_base_t> fork() const {
      return boost::shared_ptr<f_calc_function_base_t>(
        new f_calc_function_ed(*this));
    }

    virtual FloatType get_observable() const {
      return observables[index];
    }
    virtual std::complex<FloatType> get_f_calc() const {
      return f_calc[index];
    }
    virtual af::const_ref<FloatType> get_grad_observable() const {
      typedef af::versa_plain<FloatType> one_dim_type;
      typedef typename one_dim_type::accessor_type one_dim_accessor_type;
      one_dim_accessor_type a(design_matrix.accessor().n_columns());
      return af::const_ref<FloatType>(&design_matrix(index, 0), a);
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
