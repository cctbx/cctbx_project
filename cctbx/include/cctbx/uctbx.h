/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Refactored (R.W. Grosse-Kunstleve)
     2001 Jul: Merged from CVS branch sgtbx_special_pos (rwgk)
     2001 Apr: SourceForge release (R.W. Grosse-Kunstleve)
 */

/*! \file
    Toolbox for the handling of unit cell parameters.
 */

#ifndef CCTBX_UCTBX_H
#define CCTBX_UCTBX_H

#include <cmath>
#include <scitbx/constants.h>
#include <scitbx/sym_mat3.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/small.h>
#include <scitbx/array_family/shared.h>
#include <cctbx/math/utils.h>
#include <cctbx/coordinates.h>
#include <cctbx/miller.h>

namespace cctbx {

  // forward declaration
  namespace sgtbx { class rot_mx; }

  //! Shorthand for default vec3 type in unit cell toolbox.
  typedef scitbx::vec3<double> uc_vec3;
  //! Shorthand for default mat3 type in unit cell toolbox.
  typedef scitbx::mat3<double> uc_mat3;
  //! Shorthand for default sym_mat3 type in unit cell toolbox.
  typedef scitbx::sym_mat3<double> uc_sym_mat3;

  //! Unit Cell Toolbox namespace.
  namespace uctbx {

  //! Conversion of d-spacing measures.
  inline double d_star_sq_as_stol_sq(double d_star_sq)
  {
    return d_star_sq * .25;
  }

  //! Conversion of d-spacing measures.
  inline double d_star_sq_as_two_stol(double d_star_sq)
  {
    return std::sqrt(d_star_sq);
  }

  //! Conversion of d-spacing measures.
  inline double d_star_sq_as_stol(double d_star_sq)
  {
    return std::sqrt(d_star_sq) * .5;
  }

  //! Conversion of d-spacing measures.
  inline double d_star_sq_as_d(double d_star_sq)
  {
    if (d_star_sq == 0.) return -1.;
    return 1. / std::sqrt(d_star_sq);
  }

  //! Conversion of d-spacing measures.
  inline double d_star_sq_as_two_theta(double d_star_sq, double wavelength,
                                       bool deg = false)
  {
    double result = 2. * std::asin(d_star_sq_as_stol(d_star_sq) * wavelength);
    if (deg) return scitbx::rad_as_deg(result);
    return result;
  }

  //! Class for the handling of unit cell information.
  /*! All angles are in degrees.
      <p>
      The PDB convention for orthogonalization and fractionalization
      of coordinates is used:
      <pre>
      Crystallographic Basis: D = {a,b,c}
      Cartesian Basis:        C = {i,j,k}
      i || a
      j is in (a,b) plane
      k = i x j
      </pre>
    */
  class unit_cell
  {
    public:
      //! Default constructor. Some data members are not initialized!
      /*! volume() of default constructed instances == 0.
          This may be used to test if a unit_cell instance is valid.
          <p>
          Not available in Python.
       */
      unit_cell() : volume_(0) {}

      //! Constructor using parameters (a, b, c, alpha, beta, gamma).
      /*! Missing lengths are set to 1, missing angles to 90.
          <p>
          If is_metrical_matrix == true, the input is assumed to
          consist of the coefficients of the matrical matrix, of
          which there must be exactly six in the order:
          a*a, b*b, c*c, a*b*cos(gamma), a*c*cos(beta), b*c*cos(alpha)
          <p>
          See also: paramters(), metrical_matrix()
       */
      explicit
      unit_cell(af::small<double, 6> const& parameters,
                bool is_metrical_matrix = false);

      //! Constructor using parameters (a, b, c, alpha, beta, gamma).
      explicit
      unit_cell(af::double6 const& parameters);

      //! Constructor using parameters derived from a metrical matrix.
      /*! The metrical matrix is defined as:
         <pre>
         ( a*a,            a*b*cos(gamma), a*c*cos(beta)  )
         ( a*b*cos(gamma), b*b,            b*c*cos(alpha) )
         ( a*c*cos(beta),  b*c*cos(alpha), c*c            )
         </pre>
       */
      explicit
      unit_cell(uc_sym_mat3 const& metrical_matrix);

      //! Access to parameters.
      af::double6 const&
      parameters() const { return params_; }

      //! Access to reciprocal cell parameters.
      af::double6 const&
      reciprocal_parameters() const { return r_params_; }

      //! Access to metrical matrix.
      uc_sym_mat3 const&
      metrical_matrix() const { return metr_mx_; }

      //! Access to reciprocal cell metrical matrix.
      uc_sym_mat3 const&
      reciprocal_metrical_matrix() const { return r_metr_mx_; }

      //! Volume of the unit cell.
      double
      volume() const { return volume_; }

      //! Corresponding reciprocal cell.
      unit_cell
      reciprocal() const;

      //! Length^2 of the longest lattice vector in the unit cell.
      double
      longest_vector_sq() const;

      //! Comparison of unit cell parameters.
      bool
      is_similar_to(unit_cell const& other,
                    double relative_length_tolerance=0.01,
                    double absolute_angle_tolerance=1.) const;

      //! Matrix for the conversion of cartesian to fractional coordinates.
      /*! x(fractional) = matrix * x(cartesian). */
      uc_mat3 const& fractionalization_matrix() const { return frac_; }

      //! Matrix for the conversion of fractional to cartesian coordinates.
      /*! x(cartesian) = matrix * x(fractional). */
      uc_mat3 const& orthogonalization_matrix() const { return orth_; }

      //! Conversion of cartesian coordinates to fractional coordinates.
      template <class FloatType>
      fractional<FloatType>
      fractionalize(cartesian<FloatType> const& xc) const
      {
        return frac_ * xc;
      }

      //! Conversion of fractional coordinates to cartesian coordinates.
      template <class FloatType>
      cartesian<FloatType>
      orthogonalize(fractional<FloatType> const& xf) const
      {
        return orth_ * xf;
      }

      //! Length^2 of a vector of fractional coordinates.
      /*! Not available in Python.
       */
      template <class FloatType>
      FloatType
      length_sq(fractional<FloatType> const& xf) const
      {
        return orthogonalize(xf).length_sq();
      }

      //! Length of a vector of fractional coordinates.
      template <class FloatType>
      FloatType
      length(fractional<FloatType> const& xf) const
      {
        return std::sqrt(length_sq(xf));
      }

      //! Distance^2 between two vectors of fractional coordinates.
      /*! Not available in Python.
       */
      template <class FloatType>
      FloatType
      distance_sq(fractional<FloatType> const& xf,
                  fractional<FloatType> const& yf) const
      {
        return length_sq(fractional<FloatType>(xf - yf));
      }

      //! Distance between two vectors of fractional coordinates.
      template <class FloatType>
      FloatType
      distance(fractional<FloatType> const& xf,
               fractional<FloatType> const& yf) const
      {
        return length(fractional<FloatType>(xf - yf));
      }

      /*! \brief Shortest length^2 of a vector of fractional coordinates
          under application of periodicity. */
      /*! Not available in Python.
       */
      template <class FloatType>
      FloatType
      mod_short_length_sq(fractional<FloatType> const& xf) const
      {
        return length_sq(xf.mod_short());
      }

      /*! \brief Shortest length of a vector of fractional coordinates
          under application of periodicity. */
      template <class FloatType>
      FloatType
      mod_short_length(fractional<FloatType> const& xf) const
      {
        return std::sqrt(mod_short_length_sq(xf));
      }

      /*! \brief Shortest distance^2 between two vectors of fractional
          coordinates under application of periodicity. */
      /*! Not available in Python.
       */
      template <class FloatType>
      FloatType
      mod_short_distance_sq(fractional<FloatType> const& xf,
                            fractional<FloatType> const& yf) const
      {
        return mod_short_length_sq(fractional<FloatType>(xf - yf));
      }

      /*! \brief Shortest distance between two vectors of fractional
          coordinates under application of periodicity. */
      template <class FloatType>
      FloatType
      mod_short_distance(fractional<FloatType> const& xf,
                         fractional<FloatType> const& yf) const
      {
        return std::sqrt(mod_short_distance_sq(xf, yf));
      }

      /*! \brief Shortest distance^2 between all sites in xf and yf
          under application of periodicity.
       */
      /*! Not available in Python.
       */
      template <class FloatType>
      FloatType
      min_mod_short_distance_sq(
        af::const_ref<scitbx::vec3<FloatType> > const& xf,
        fractional<FloatType> const& yf) const
      {
        FloatType
          result = mod_short_distance_sq(fractional<FloatType>(xf[0]), yf);
        for(std::size_t i=1;i<xf.size();i++) {
          math::update_min(
            result, mod_short_distance_sq(fractional<FloatType>(xf[i]), yf));
        }
        return result;
      }

      /*! \brief Shortest distance between all sites in xf and yf
          under application of periodicity.
       */
      template <class FloatType>
      FloatType
      min_mod_short_distance(
        af::const_ref<scitbx::vec3<FloatType> > const& xf,
        fractional<FloatType> const& yf) const
      {
        return std::sqrt(min_mod_short_distance_sq(xf, yf));
      }

      //! Transformation (change-of-basis) of unit cell parameters.
      /*! r is the inverse of the 3x3 change-of-basis matrix
          that transforms coordinates in the old basis system to
          coodinates in the new basis system.
       */
      unit_cell
      change_basis(uc_mat3 const& r, double r_den=1.) const;

      //! Transformation (change-of-basis) of unit cell parameters.
      /*! r is the inverse of the 3x3 change-of-basis matrix
          that transforms coordinates in the old basis system to
          coodinates in the new basis system.
          <p>
          See also: sgtbx::change_of_basis::apply()
          <p>
          Not available in Python.
       */
      unit_cell
      change_basis(sgtbx::rot_mx const& c_inv_r) const;

      /*! \brief Computation of the maximum Miller indices for a given
          minimum d-spacing.
       */
      /*! d_min is the minimum d-spacing. tolerance compensates for
          rounding errors.
       */
      miller::index<>
      max_miller_indices(double d_min, double tolerance=1.e-4) const;

      //! d-spacing measure d_star_sq = 1/d^2 = (2*sin(theta)/lambda)^2.
      template <typename NumType>
      double
      d_star_sq(miller::index<NumType> const& h) const
      {
        return
            (h[0] * h[0]) * r_metr_mx_[0]
          + (h[1] * h[1]) * r_metr_mx_[1]
          + (h[2] * h[2]) * r_metr_mx_[2]
          + (2 * h[0] * h[1]) * r_metr_mx_[3]
          + (2 * h[0] * h[2]) * r_metr_mx_[4]
          + (2 * h[1] * h[2]) * r_metr_mx_[5];
      }

      //! d-spacing measure d_star_sq = 1/d^2 = (2*sin(theta)/lambda)^2.
      template <typename NumType>
      af::shared<double>
      d_star_sq(af::const_ref<miller::index<NumType> > const& h) const
      {
        af::shared<double> result(h.size(), af::init_functor_null<double>());
        for(std::size_t i=0;i<result.size();i++) {
          result[i] = d_star_sq(h[i]);
        }
        return result;
      }

      //! Maximum d_star_sq for given list of Miller indices.
      template <typename NumType>
      double
      max_d_star_sq(af::const_ref<miller::index<NumType> > const& h) const
      {
        double result = 0;
        for(std::size_t i=0;i<h.size();i++) {
          math::update_max(result, d_star_sq(h[i]));
        }
        return result;
      }

      //! Minimum and maximum d_star_sq for given list of Miller indices.
      template <typename NumType>
      af::double2
      min_max_d_star_sq(af::const_ref<miller::index<NumType> > const& h) const
      {
        af::double2 result(0, 0);
        if (h.size()) {
          result.fill(d_star_sq(h[0]));
          for(std::size_t i=1;i<h.size();i++) {
            double q = d_star_sq(h[i]);
            math::update_min(result[0], q);
            math::update_max(result[1], q);
          }
        }
        return result;
      }

      //! d-spacing measure (sin(theta)/lambda)^2 = d_star_sq/4.
      template <typename NumType>
      double
      stol_sq(miller::index<NumType> const& h) const
      {
        return d_star_sq_as_stol_sq(d_star_sq(h));
      }

      //! d-spacing measure (sin(theta)/lambda)^2 = d_star_sq/4.
      template <typename NumType>
      af::shared<double>
      stol_sq(af::const_ref<miller::index<NumType> > const& h) const
      {
        af::shared<double> result(h.size(), af::init_functor_null<double>());
        for(std::size_t i=0;i<result.size();i++) {
          result[i] = stol_sq(h[i]);
        }
        return result;
      }

      //! d-spacing measure 2*sin(theta)/lambda = 1/d = sqrt(d_star_sq).
      template <typename NumType>
      double
      two_stol(miller::index<NumType> const& h) const
      {
        return d_star_sq_as_two_stol(d_star_sq(h));
      }

      //! d-spacing measure 2*sin(theta)/lambda = 1/d = sqrt(d_star_sq).
      template <typename NumType>
      af::shared<double>
      two_stol(af::const_ref<miller::index<NumType> > const& h) const
      {
        af::shared<double> result(h.size(), af::init_functor_null<double>());
        for(std::size_t i=0;i<result.size();i++) {
          result[i] = two_stol(h[i]);
        }
        return result;
      }

      //! d-spacing measure sin(theta)/lambda = 1/(2*d) = sqrt(d_star_sq)/2.
      template <typename NumType>
      double
      stol(miller::index<NumType> const& h) const
      {
        return d_star_sq_as_stol(d_star_sq(h));
      }

      //! d-spacing measure sin(theta)/lambda = 1/(2*d) = sqrt(d_star_sq)/2.
      template <typename NumType>
      af::shared<double>
      stol(af::const_ref<miller::index<NumType> > const& h) const
      {
        af::shared<double> result(h.size(), af::init_functor_null<double>());
        for(std::size_t i=0;i<result.size();i++) {
          result[i] = stol(h[i]);
        }
        return result;
      }

      //! d-spacing measure d = 1/(2*sin(theta)/lambda).
      template <typename NumType>
      double
      d(miller::index<NumType> const& h) const
      {
        return d_star_sq_as_d(d_star_sq(h));
      }

      //! d-spacing measure d = 1/(2*sin(theta)/lambda).
      template <typename NumType>
      af::shared<double>
      d(af::const_ref<miller::index<NumType> > const& h) const
      {
        af::shared<double> result(h.size(), af::init_functor_null<double>());
        for(std::size_t i=0;i<result.size();i++) {
          result[i] = d(h[i]);
        }
        return result;
      }

      //! Diffraction angle 2-theta, given wavelength.
      template <typename NumType>
      double
      two_theta(
        miller::index<NumType> const& h,
        double wavelength,
        bool deg = false) const
      {
        return d_star_sq_as_two_theta(d_star_sq(h), wavelength, deg);
      }

      //! Diffraction angle 2-theta, given wavelength.
      template <typename NumType>
      af::shared<double>
      two_theta(
        af::const_ref<miller::index<NumType> > const& h,
        double wavelength,
        bool deg = false) const
      {
        af::shared<double> result(h.size(), af::init_functor_null<double>());
        for(std::size_t i=0;i<result.size();i++) {
          result[i] = two_theta(h[i], wavelength, deg);
        }
        return result;
      }

    protected:
      void init_volume();
      void init_reciprocal();
      void init_metrical_matrices();
      void init_orth_and_frac_matrices();
      void initialize();

      af::double6 params_;
      af::double3 sin_ang_;
      af::double3 cos_ang_;
      double volume_;
      uc_sym_mat3 metr_mx_;
      af::double6 r_params_;
      af::double3 r_sin_ang_;
      af::double3 r_cos_ang_;
      uc_sym_mat3 r_metr_mx_;
      uc_mat3 frac_;
      uc_mat3 orth_;
      mutable double longest_vector_sq_;

      // used by reciprocal()
      unit_cell(
        af::double6 const& params,
        af::double3 const& sin_ang,
        af::double3 const& cos_ang,
        double volume,
        uc_sym_mat3 const& metr_mx,
        af::double6 const& r_params,
        af::double3 const& r_sin_ang,
        af::double3 const& r_cos_ang,
        uc_sym_mat3 const& r_metr_mx);
  };

  /*! \brief Helper class for optimizing d_star_sq computations in
      loops over a grid of Miller indices.
   */
  template <typename FloatType>
  class incremental_d_star_sq
  {
    public:
      //! Default contructor. Some data members are not initialized!
      incremental_d_star_sq() {}

      //! Initialization from unit_cell object.
      /*! This copies the elements of the metrical matrix.
       */
      incremental_d_star_sq(unit_cell const& ucell)
      {
        initialize(ucell.reciprocal_metrical_matrix());
      }

      //! Stores h0 and performs computations that only involve h0.
      void update0(int h0)
      {
        h0_ = h0;
        im0_ = (h0_ * h0_) * r_g00_;
      }

      //! Stores h1 and performs computations that only involve h0 and h1.
      void update1(int h1)
      {
        h1_ = h1;
        im1_ = im0_ + (h1_ * h1_) * r_g11_
                    + (2 * h0_ * h1_) * r_g01_;
      }

      //! Returns d_star_sq using (h0,h1,h2).
      FloatType get(int h2)
      {
        return im1_ + (h2 * h2) * r_g22_
                    + (2 * h0_ * h2) * r_g02_
                    + (2 * h1_ * h2) * r_g12_;
      }

    protected:
      FloatType r_g00_, r_g11_, r_g22_, r_g01_, r_g02_, r_g12_;
      int h0_, h1_;
      FloatType im0_, im1_;

      void initialize(uc_sym_mat3 const& r_g)
      {
        r_g00_ = r_g[0];
        r_g11_ = r_g[1];
        r_g22_ = r_g[2];
        r_g01_ = r_g[3];
        r_g02_ = r_g[4];
        r_g12_ = r_g[5];
      }
  };

}} // namespace cctbx::uctbx

#endif // CCTBX_UCTBX_H
