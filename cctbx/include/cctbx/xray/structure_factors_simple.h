#ifndef CCTBX_XRAY_STRUCTURE_FACTORS_SIMPLE_H
#define CCTBX_XRAY_STRUCTURE_FACTORS_SIMPLE_H

/* Reference implementation of structure factor calculation.
   For regression tests only.
   Simple code but slow.
   Do not use for applications.
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include <cctbx/xray/scattering_dictionary.h>
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
        if (scatterer.anisotropic_flag) {
          dw = adptbx::debye_waller_factor_u_star(hr, scatterer.u_star);
        }
        else {
          dw = adptbx::debye_waller_factor_u_iso(d_star_sq/4, scatterer.u_iso);
        }
        f_calc += e_j_phase * dw;
      }
      f_calc *= (f0 + c_t(scatterer.fp, scatterer.fdp)) * scatterer.weight();
    }

    std::complex<float_type> f_calc;
  };

  template <typename ScattererType=scatterer<> >
  struct simple_one_h
  {
    typedef typename ScattererType::float_type float_type;

    simple_one_h(
      uctbx::unit_cell const& unit_cell,
      sgtbx::space_group const& space_group,
      miller::index<> const& h,
      af::const_ref<ScattererType> const& scatterers,
      scattering_dictionary const& scattering_dict)
    :
      f_calc(0,0)
    {
      typedef scattering_dictionary::dict_type dict_type;
      typedef dict_type::const_iterator dict_iter;
      typedef float_type f_t;
      f_t d_star_sq = unit_cell.d_star_sq(h);
      dict_type const& scd = scattering_dict.dict();
      for(dict_iter di=scd.begin();di!=scd.end();di++) {
        f_t f0 = di->second.gaussian.at_d_star_sq(d_star_sq);
        af::const_ref<std::size_t>
          member_indices = di->second.member_indices.const_ref();
        for(std::size_t mi=0;mi<member_indices.size();mi++) {
          ScattererType const& scatterer = scatterers[member_indices[mi]];
          //CCTBX_ASSERT(scatterer.scatting_type == di->first);
          simple_one_h_one_scatterer<ScattererType> sf(
            space_group,
            h,
            d_star_sq,
            scatterer,
            f0);
          f_calc += sf.f_calc;
        }
      }
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
        scattering_dictionary const& scattering_dict)
      {
        CCTBX_ASSERT(scattering_dict.n_scatterers() == scatterers.size());
        f_calc_.reserve(miller_indices.size());
        for(std::size_t i=0;i<miller_indices.size();i++) {
          f_calc_.push_back(simple_one_h<ScattererType>(
            unit_cell,
            space_group,
            miller_indices[i],
            scatterers,
            scattering_dict).f_calc);
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
