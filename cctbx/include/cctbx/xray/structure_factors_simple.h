#ifndef CCTBX_XRAY_STRUCTURE_FACTORS_SIMPLE_H
#define CCTBX_XRAY_STRUCTURE_FACTORS_SIMPLE_H

/* Reference implementation of structure factor calculation.
   For regression tests only.
   Simple code but slow.
   Do not use for applications.
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include <cctbx/xray/scatterer.h>
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
      ScattererType const& scatterer)
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
        if (scatterer.anisotropic_flag) {
          dw = adptbx::debye_waller_factor_u_star(hr, scatterer.u_star);
        }
        else {
          dw = adptbx::debye_waller_factor_u_iso(d_star_sq/4, scatterer.u_iso);
        }
        f_calc += e_j_phase * dw;
      }
      f_t f0 = scatterer.caasf.at_d_star_sq(d_star_sq);
      f_calc *= (f0 + c_t(scatterer.fp, scatterer.fdp)) * scatterer.weight();
    }

    std::complex<float_type> f_calc;
  };

  template <typename ScattererType=scatterer<> >
  struct simple_one_scatterer
  {
    typedef typename ScattererType::float_type float_type;

    simple_one_scatterer(
      uctbx::unit_cell const& unit_cell,
      sgtbx::space_group const& space_group,
      af::const_ref<miller::index<> > const& miller_indices,
      ScattererType const& scatterer,
      af::ref<std::complex<float_type> > const& f_calc)
    {
      typedef float_type f_t;
      for(std::size_t i=0;i<miller_indices.size();i++) {
        miller::index<> const& h = miller_indices[i];
        f_t d_star_sq = unit_cell.d_star_sq(h);
        simple_one_h_one_scatterer<ScattererType> sf(
          space_group,
          h,
          d_star_sq,
          scatterer);
        f_calc[i] += sf.f_calc;
      }
    }
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
        af::const_ref<ScattererType> const& scatterers)
      {
        f_calc_.resize(miller_indices.size());
        for(std::size_t i=0;i<scatterers.size();i++) {
          ScattererType const& scatterer = scatterers[i];
          simple_one_scatterer<ScattererType> sf(
            unit_cell,
            space_group,
            miller_indices,
            scatterer,
            f_calc_.ref());
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
