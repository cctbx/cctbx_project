#ifndef CCTBX_XRAY_SCATTERER_H
#define CCTBX_XRAY_SCATTERER_H

#include <cctbx/sgtbx/site_symmetry.h>
#include <cctbx/adptbx.h>

namespace cctbx {
  //! X-ray scatterer and structure factor namespace.
namespace xray {

  /*! \brief Information about an atom that is needed for a
      structure factor calculation.
   */
  /*! Constructors are provided for scatterers with isotropic
      and anisotropic displacement parameters (temperature factors).
      <p>
      The anisotropic displacement parameters have to be
      "u_star." Conversions between different conventions for
      the representation of anisotropic displacement
      parameters are provided by cctbx::adptbx.
      <p>
      The apply_symmetry() member function has to be called before
      the scatterer can be used in a structure factor calculation.
   */
  template <typename FloatType=double,
            typename LabelType=std::string,
            typename ScatteringTypeType=std::string>
  class scatterer
  {
    public:
      //! Facilitates meta-programming.
      typedef FloatType float_type;
      //! Facilitates meta-programming.
      typedef LabelType label_type;
      //! Facilitates meta-programming.
      typedef ScatteringTypeType scattering_type_type;

      //! Default constructor. Data members are not initialized!
      scatterer() {}

      //! Initialization with isotropic displacement parameter.
      scatterer(LabelType const& label_,
                fractional<FloatType> const& site_,
                FloatType const& u_iso_,
                FloatType const& occupancy_,
                ScatteringTypeType const& scattering_type_,
                FloatType fp_,
                FloatType fdp_)
      :
        label(label_),
        scattering_type(scattering_type_),
        fp(fp_),
        fdp(fdp_),
        site(site_),
        occupancy(occupancy_),
        anisotropic_flag(false),
        u_iso(u_iso_),
        u_star(-1,-1,-1,-1,-1,-1),
        multiplicity_(0),
        weight_without_occupancy_(0)
      {}

      //! Initialization with anisotropic displacement parameters.
      scatterer(LabelType const& label_,
                fractional<FloatType> const& site_,
                scitbx::sym_mat3<FloatType> const& u_star_,
                FloatType const& occupancy_,
                ScatteringTypeType const& scattering_type_,
                FloatType fp_,
                FloatType fdp_)
      :
        label(label_),
        scattering_type(scattering_type_),
        fp(fp_),
        fdp(fdp_),
        site(site_),
        occupancy(occupancy_),
        anisotropic_flag(true),
        u_iso(-1),
        u_star(u_star_),
        multiplicity_(0),
        weight_without_occupancy_(0)
      {}

      //! Direct access to label.
      LabelType label;

      //! Direct access to the scattering type.
      /*! See also: eltbx::xray_scattering::it1992,
                    eltbx::xray_scattering::wk1995
       */
      ScatteringTypeType scattering_type;

      //! Direct access to f-prime.
      /*! f-prime is the dispersive contribution to the scattering
          factor.
       */
      FloatType fp;

      //! Direct access to f-double-prime.
      /*! f-double-prime is the anomalous contribution to the scattering
          factor.
       */
      FloatType fdp;

      //! Direct access to fractional coordinates.
      /*! See also: apply_symmetry(), apply_symmetry_site()
       */
      fractional<FloatType> site;

      //! Direct access to occupancy factor.
      FloatType occupancy;

      //! Direct access to flag indicating anisotropic displacement parameters.
      bool anisotropic_flag;

      //! Direct access to isotropic displacement parameter.
      /*! Defined only if anisotropic_flag == false.
          <p>
          Conversions between isotropic and anisotropic displacement
          parameters are provided by cctbx::adptbx.
       */
      FloatType u_iso;

      //! Direct access to anisotropic displacement parameters.
      /*! Defined only if anisotropic_flag == true.
          <p>
          Conversions between isotropic and anisotropic displacement
          parameters are provided by cctbx::adptbx.
          <p>
          See also: apply_symmetry(), apply_symmetry_u_star()
       */
      scitbx::sym_mat3<FloatType> u_star;

      /*! \brief Computes multiplicity(), weight_without_occupancy(),
          weight() and symmetry-averaged anisotropic displacement parameters.
       */
      /*! This member function must be called before the
          scatterer is used in a structure factor calculation.

          See also:
            apply_symmetry_site,
            apply_symmetry_u_star,
            cctbx::sgtbx::site_symmetry,
            cctbx::sgtbx::site_symmetry::exact_site,
            cctbx::sgtbx::site_symmetry_ops::average_u_star
       */
      sgtbx::site_symmetry
      apply_symmetry(
        uctbx::unit_cell const& unit_cell,
        sgtbx::space_group const& space_group,
        double min_distance_sym_equiv=0.5,
        double u_star_tolerance=0,
        bool assert_is_positive_definite=false,
        bool assert_min_distance_sym_equiv=true);

      //! Apply previously determined site symmetry to site.
      /*! Shorthand for: site = site_symmetry_ops.special_op() * site
       */
      void
      apply_symmetry_site(sgtbx::site_symmetry_ops const& site_symmetry_ops)
      {
        site = site_symmetry_ops.special_op() * site;
      }

      //! Apply previously determined site symmetry to u_star.
      /*! For scatterers with anisotropic displacement
          parameters, the symmetry-averaged u_star is determined
          using cctbx::sgtbx::site_symmetry_ops::average_u_star .
          If u_star_tolerance is greater than zero an
          exception is thrown if the discrepancy between
          components of u_star before and after the
          application of the site symmetry is greater than
          u_star_tolerance.

          If assert_is_positive_definite == true,
          for scatterers with anisotropic displacement
          parameters it is tested if the symmetry-averaged
          u_star tensor is positive definite. An exception
          is thrown if this is not the case.
          adptbx::eigenvalue_filtering is used to reset
          negative eigenvalues of u_star to zero.
          cctbx::sgtbx::site_symmetry::average_u_star
          is applied again to compute the final u_star.

          This function has no effect if anisotropic_flag == false.
       */
      void
      apply_symmetry_u_star(
        uctbx::unit_cell const& unit_cell,
        sgtbx::site_symmetry_ops const& site_symmetry_ops,
        double u_star_tolerance=0,
        bool assert_is_positive_definite=false,
        bool assert_min_distance_sym_equiv=true);

      //! Access to multiplicity computed by apply_symmetry().
      int
      multiplicity() const { return multiplicity_; }

      //! Access to "weight_without_occupancy" computed by apply_symmetry().
      /*! weight_without_occupancy is defined as
              multiplicity() / space_group.order_z(),
          with space_group as passed to apply_symmetry().
          weight_without_occupancy() is used in the computation
          of structure factor derivatives.
       */
      FloatType
      weight_without_occupancy() const { return weight_without_occupancy_; }

      //! Access to "weight" computed by apply_symmetry().
      /*! The weight is defined as
              occupancy * multiplicity() / space_group.order_z(),
          with space_group as passed to apply_symmetry().
          The weight() is used in structure factor calculations.
       */
      FloatType
      weight() const { return weight_without_occupancy_ * occupancy; }

      //! Helper function for object serialization (Python pickle).
      /*! For internal use only.
       */
      void
      setstate(int multiplicity,
               FloatType weight_without_occupancy)
      {
        multiplicity_ = multiplicity;
        weight_without_occupancy_ = weight_without_occupancy;
      }

    protected:
      int multiplicity_;
      FloatType weight_without_occupancy_;
  };

  template <typename FloatType,
            typename LabelType,
            typename ScatteringTypeType>
  sgtbx::site_symmetry
  scatterer<FloatType, LabelType, ScatteringTypeType>
  ::apply_symmetry(
    uctbx::unit_cell const& unit_cell,
    sgtbx::space_group const& space_group,
    double min_distance_sym_equiv,
    double u_star_tolerance,
    bool assert_is_positive_definite,
    bool assert_min_distance_sym_equiv)
  {
    sgtbx::site_symmetry site_symmetry(
      unit_cell,
      space_group,
      site,
      min_distance_sym_equiv,
      assert_min_distance_sym_equiv);
    multiplicity_ = site_symmetry.multiplicity();
    weight_without_occupancy_ = FloatType(multiplicity_)
                              / space_group.order_z();
    site = site_symmetry.exact_site();
    apply_symmetry_u_star(
      unit_cell,
      site_symmetry,
      u_star_tolerance,
      assert_is_positive_definite,
      assert_min_distance_sym_equiv);
    return site_symmetry;
  }

  template <typename FloatType,
            typename LabelType,
            typename ScatteringTypeType>
  void
  scatterer<FloatType, LabelType, ScatteringTypeType>
  ::apply_symmetry_u_star(
    uctbx::unit_cell const& unit_cell,
    sgtbx::site_symmetry_ops const& site_symmetry_ops,
    double u_star_tolerance,
    bool assert_is_positive_definite,
    bool assert_min_distance_sym_equiv)
  {
    if (anisotropic_flag) {
      if (u_star_tolerance > 0.) {
        CCTBX_ASSERT(
          site_symmetry_ops.is_compatible_u_star(u_star, u_star_tolerance));
      }
      u_star = site_symmetry_ops.average_u_star(u_star);
      scitbx::sym_mat3<FloatType>
        u_cart = adptbx::u_star_as_u_cart(unit_cell, u_star);
      if (assert_is_positive_definite) {
        CCTBX_ASSERT(adptbx::is_positive_definite(u_cart));
      }
      u_cart = adptbx::eigenvalue_filtering(u_cart);
      u_star = adptbx::u_cart_as_u_star(unit_cell, u_cart);
      u_star = site_symmetry_ops.average_u_star(u_star);
    }
  }

}} // namespace cctbx::xray

#endif // CCTBX_XRAY_SCATTERER_H
