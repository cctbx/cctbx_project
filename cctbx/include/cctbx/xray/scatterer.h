/* Copyright (c) 200f_calc-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Refactored (rwgk)
     2001 Oct: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_XRAY_SCATTERER_H
#define CCTBX_XRAY_SCATTERER_H

#include <cctbx/eltbx/caasf.h>
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
  template <typename FloatType = double,
            typename CaasfType = eltbx::caasf::wk1995,
            typename LabelType = std::string>
  class scatterer
  {
    public:
      //! Facilitates meta-programming.
      typedef FloatType float_type;
      //! Facilitates meta-programming.
      typedef CaasfType caasf_type;
      //! Facilitates meta-programming.
      typedef LabelType label_type;

      //! Default constructor. Data members are not initialized!
      scatterer() {}

      //! Initialization with isotropic displacement parameter.
      scatterer(LabelType const& label_,
                fractional<FloatType> const& site_,
                FloatType const& u_iso_,
                FloatType const& occupancy_,
                CaasfType const& caasf_,
                std::complex<FloatType> const& fp_fdp_)
      :
        label(label_),
        caasf(caasf_),
        fp_fdp(fp_fdp_),
        site(site_),
        occupancy(occupancy_),
        anisotropic_flag(false),
        u_iso(u_iso_),
        u_star(-1,-1,-1,-1,-1,-1),
        multiplicity_(0),
        weight_(0)
      {}

      //! Initialization with anisotropic displacement parameters.
      scatterer(LabelType const& label_,
                fractional<FloatType> const& site_,
                scitbx::sym_mat3<FloatType> const& u_star_,
                FloatType const& occupancy_,
                CaasfType const& caasf_,
                std::complex<FloatType> const& fp_fdp_)
      :
        label(label_),
        caasf(caasf_),
        fp_fdp(fp_fdp_),
        site(site_),
        occupancy(occupancy_),
        anisotropic_flag(true),
        u_iso(-1),
        u_star(u_star_),
        multiplicity_(0),
        weight_(0)
      {}

      //! Direct access to label.
      LabelType label;

      //! Direct access to analytical approximation to the scattering factor.
      /*! See also: eltbx::caasf::it1992,
                    eltbx::caasf::wk1995
       */
      CaasfType caasf;

      //! Direct access to f-prime and f-double-prime.
      /*! f-prime is the dispersive contribution to the scattering
          factor. f-double-prime is the anomalous contribution.
       */
      std::complex<FloatType> fp_fdp;

      //! Direct access to fractional coordinates.
      /*! See also: apply_symmetry()
       */
      fractional<FloatType> site;

      //! Direct access to occupancy factor.
      /*! See also: apply_symmetry()
       */
      FloatType occupancy;

      //! Direct access to flag indicating anisotropic displacement parameters.
      bool anisotropic_flag;

      //! Direct access to isotropic displacement parameter.
      /*! Defined only if anisotropic_flag == false.
          <p>
          Conversions between isotropic and anisotropic displacement
          parameters are provided by cctbx::adptbx.
          <p>
          See also: apply_symmetry()
       */
      FloatType u_iso;

      //! Direct access to anisotropic displacement parameters.
      /*! Defined only if anisotropic_flag == true.
          <p>
          Conversions between isotropic and anisotropic displacement
          parameters are provided by cctbx::adptbx.
          <p>
          See also: apply_symmetry()
       */
      scitbx::sym_mat3<FloatType> u_star;

      /*! \brief Computes multiplicity(), weight() and symmetry-averaged
           anisotropic displacement parameters.
       */
      /*! This member function must be called before the
          scatterer is used in a structure factor calculation
          and after the last change of site and u_star.
          <p>
          See also: cctbx::sgtbx::site_symmetry,
                    cctbx::sgtbx::site_symmetry::average_u_star
          <p>
          For scatterers with anisotropic displacement
          parameters, the symmetry-averaged u_star is determined
          using cctbx::sgtbx::site_symmetry::average_u_star .
          If u_star_tolerance is greater than zero an
          exception is thrown if the discrepancy between
          components of u_star before and after the
          application of the site symmetry is greater than
          u_star_tolerance.
          <p>
          If assert_is_positive_definite == true,
          for scatterers with anisotropic displacement
          parameters it is tested if the symmetry-averaged
          u_star tensor is positive definite. An exception
          is thrown if this is not the case.
       */
      sgtbx::site_symmetry
      apply_symmetry(uctbx::unit_cell const& unit_cell,
                     sgtbx::space_group const& space_group,
                     double min_distance_sym_equiv=0.5,
                     double u_star_tolerance=0.1,
                     bool assert_is_positive_definite=true,
                     bool assert_min_distance_sym_equiv=true);

      //! Updates the weight().
      /*! This member function must be called if the occupancy
          factor is changed after calling apply_symmetry(),
          and calling apply_symmetry() again is not desired
          for some reason.
       */
      void
      update_weight(sgtbx::space_group const& space_group);

      //! Access to multiplicity computed by apply_symmetry().
      int
      multiplicity() const { return multiplicity_; }

      //! Access to "weight" computed by apply_symmetry().
      /*! The weight is defined as
              occupancy() * multiplicity() / space_group.order_z(),
          with space_group as passed to apply_symmetry().
          The weight is used in structure factor calculations.
       */
      FloatType
      weight() const { return weight_; }

      //! Helper function for object serialization (Python pickle).
      /*! For internal use only.
       */
      void
      setstate(std::string const& caasf_label,
               int multiplicity,
               FloatType weight)
      {
        caasf = caasf_type(caasf_label);
        multiplicity_ = multiplicity;
        weight_ = weight;
      }

    protected:
      int multiplicity_;
      FloatType weight_;
  };

  template <typename FloatType,
            typename CaasfType,
            typename LabelType>
  sgtbx::site_symmetry
  scatterer<FloatType, CaasfType, LabelType>
  ::apply_symmetry(uctbx::unit_cell const& unit_cell,
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
    site = site_symmetry.exact_site();
    multiplicity_ = site_symmetry.multiplicity();
    weight_ = occupancy * multiplicity_ / space_group.order_z();
    if (anisotropic_flag) {
      if (u_star_tolerance > 0.) {
        CCTBX_ASSERT(
          site_symmetry.is_compatible_u_star(u_star, u_star_tolerance));
      }
      u_star = site_symmetry.average_u_star(u_star);
      if (assert_is_positive_definite) {
        CCTBX_ASSERT(adptbx::is_positive_definite(u_star));
      }
    }
    return site_symmetry;
  }

  template <typename FloatType,
            typename CaasfType,
            typename LabelType>
  void
  scatterer<FloatType, CaasfType, LabelType>
  ::update_weight(sgtbx::space_group const& space_group)
  {
    weight_ = occupancy * multiplicity_ / space_group.order_z();
  }

}} // namespace cctbx::xray

#endif // CCTBX_XRAY_SCATTERER_H
