#ifndef CCTBX_XRAY_HR_HT_CACHE_H
#define CCTBX_XRAY_HR_HT_CACHE_H

#include <cctbx/math/cos_sin_table.h>
#include <cctbx/sgtbx/miller_ops.h>

namespace cctbx { namespace xray { namespace structure_factors {

  /// A pair \f$(hR, ht)\f$ where \f$(R \mid t)\f$ is a symmetry element.
  template <typename FloatType>
  struct hr_ht_group
  {
    hr_ht_group(
      miller::index<> const& hr_,
      FloatType const& ht_)
    :
      hr(hr_),
      ht(ht_)
    {}

    /// \f$hR\f$
    miller::index<> hr;

    /// \f$ht\f$
    FloatType ht;
  };


  /** \brief The miller indices \f$h R\f$ and the scalars \f$h t\f$
       where \f$(R \mid t)\f$ spans the elements of a space group \f$G\f$.

    A centric space group is partition as the coset decomposition
    \f$G = G' \bigcup \bar{1}G''\f$ and only the elements of G' are considered.
   */
  template <typename FloatType>
  struct hr_ht_cache
  {
    typedef FloatType f_t;

    /// Construct the cache using std::cos and std::sin.
    hr_ht_cache(sgtbx::space_group const& space_group,
                miller::index<> const& h)
    {
      math::cos_sin_exact<f_t> cos_sin;
      init(cos_sin, space_group, h);
    }

    /// Construct the cache using the functor cos_sin
    /** That functor implements \f$\theta \mapsto \exp i 2\pi \theta\f$ */
    template<class CosSinType>
    hr_ht_cache(CosSinType const &cos_sin,
                sgtbx::space_group const& space_group,
                miller::index<> const& h)
    {
      init(cos_sin, space_group, h);
    }

    /// Implementation of the constructors
    template <class CosSinType>
    void init(CosSinType const &cos_sin,
              sgtbx::space_group const& space_group,
              miller::index<> const& h)
    {
      ltr_factor = space_group.n_ltr();
      is_centric = space_group.is_centric();
      f_t t_den = static_cast<f_t>(space_group.t_den());
      if (!is_centric) {
        h_inv_t = -1;
        is_origin_centric = false;
      }
      else {
        h_inv_t = static_cast<f_t>(h * space_group.inv_t()) / t_den;
        is_origin_centric = (h_inv_t == 0);
      }
      for(std::size_t i_smx=0;i_smx<space_group.n_smx();i_smx++) {
        sgtbx::rt_mx const& s = space_group.smx(i_smx);
        groups.push_back(hr_ht_group<f_t>(
          h * s.r(),
          static_cast<f_t>(h * s.t()) / t_den));
      }
      if (is_centric) {
        f_h_inv_t = !is_origin_centric ? cos_sin(h_inv_t) : 1;
      }
    }

    /// Whether the space group is centric
    bool is_centric;

    /// Whether the reflection is origin centric in the given space group
    /** i.e. whether h_inv_t is zero */
    bool is_origin_centric;

    /// Order of the centring translation subgroup
    f_t ltr_factor;

    /// \f$h t_{\bar{1}}\f$ where \f$t_{\bar{1}}\f$ is the translation part
    /// of the inversion.
    /** Only defined for a centric space-group and zero for an origin centric one
     */
    f_t h_inv_t;

    /// \f$ \exp i 2 \pi h t_{\bar{1}} \f$
    /** Only defined for a centric space group and 1 for an origin centric one */
    std::complex<f_t> f_h_inv_t;

    /// The list of pairs \f$(hR, ht)\f$.
    af::small<hr_ht_group<f_t>, sgtbx::n_max_repr_rot_mx> groups;
  };

}}} // namespace cctbx::xray::structure_factors

#endif // CCTBX_XRAY_HR_HT_CACHE_H
