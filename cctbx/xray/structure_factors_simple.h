#ifndef CCTBX_XRAY_STRUCTURE_FACTORS_SIMPLE_H
#define CCTBX_XRAY_STRUCTURE_FACTORS_SIMPLE_H

/* Reference implementation of structure factor calculation.
   For regression tests only.
   Simple code but slow.
   Do not use for applications.
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include <cctbx/xray/scattering_type_registry.h>
#include <cctbx/sgtbx/miller_ops.h>

namespace cctbx { namespace xray { namespace structure_factors {

  template <typename ScattererType=scatterer<> >
  struct simple_one_h_one_scatterer
  {
    typedef typename ScattererType::float_type float_type;

    simple_one_h_one_scatterer(
      sgtbx::space_group const& space_group,
      miller::index<> const& h,
      float_type d_star_sq,
      ScattererType const& scatterer,
      float_type f0)
    :
      f_calc(0,0)
    {
      using scitbx::constants::two_pi;
      typedef float_type f_t;
      typedef std::complex<f_t> c_t;
      f_t dw;
      for(std::size_t i_smx=0;i_smx<space_group.order_z();i_smx++) {
        sgtbx::rt_mx s = space_group(i_smx);
        miller::index<> hr = h * s.r();
        f_t hrx = hr * scatterer.site;
        f_t ht = f_t(h * s.t()) / space_group.t_den();
        f_t phase = two_pi * (hrx + ht);
        c_t e_j_phase(std::cos(phase), std::sin(phase));
        dw = 1.0;
        if (scatterer.flags.use_u_aniso()) {
          dw *= adptbx::debye_waller_factor_u_star(hr, scatterer.u_star);
          if (scatterer.anharmonic_adp) {
            e_j_phase *= scatterer.anharmonic_adp->calculate(hr);
          }
        }
        if (scatterer.flags.use_u_iso()) {
          dw *= adptbx::debye_waller_factor_u_iso(d_star_sq/4, scatterer.u_iso);
        }
        f_calc += e_j_phase * dw;
      }
      f_calc *= (f0 + c_t(scatterer.fp, scatterer.fdp)) * scatterer.weight();
    }

    std::complex<float_type> f_calc;
  };

  template <typename ScattererType=scatterer<> >
  class simple
  {
    public:
      typedef ScattererType scatterer_type;
      typedef typename ScattererType::float_type float_type;

      simple() {}

      simple(
        uctbx::unit_cell const& unit_cell,
        sgtbx::space_group const& space_group,
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<ScattererType> const& scatterers,
        xray::scattering_type_registry const& scattering_type_registry)
      {
        f_calc_.reserve(miller_indices.size());
        af::shared<std::size_t>
          scattering_type_indices_memory
            = scattering_type_registry.unique_indices(scatterers);
        af::const_ref<std::size_t>
          scattering_type_indices = scattering_type_indices_memory.const_ref();
        for(std::size_t i=0;i<miller_indices.size();i++) {
          miller::index<> h = miller_indices[i];
          double d_star_sq = unit_cell.d_star_sq(h);
          af::shared<double>
            form_factors_memory
              = scattering_type_registry.unique_form_factors_at_d_star_sq(
                  d_star_sq);
          af::const_ref<double> form_factors = form_factors_memory.const_ref();
          std::complex<float_type> f_calc(0,0);
          for(std::size_t i_sc=0;i_sc<scatterers.size();i_sc++) {
            float_type f0 = form_factors[scattering_type_indices[i_sc]];
            f_calc += simple_one_h_one_scatterer<ScattererType>(
              space_group,
              h,
              d_star_sq,
              scatterers[i_sc],
              f0).f_calc;
          }
          f_calc_.push_back(f_calc);
        }
      }

      af::shared<std::complex<float_type> >
      f_calc() const { return f_calc_; }

    protected:
      af::shared<std::complex<float_type> > f_calc_;
  };

}}} // namespace cctbx::xray::structure_factors

#endif // DOXYGEN_SHOULD_SKIP_THIS

#endif // CCTBX_XRAY_STRUCTURE_FACTORS_SIMPLE_H
