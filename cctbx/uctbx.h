#ifndef CCTBX_UCTBX_H
#define CCTBX_UCTBX_H

#include <cmath>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/optional.hpp>
#include <scitbx/constants.h>
#include <scitbx/sym_mat3.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/small.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/tiny_mat_ref.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/math/utils.h>
#include <scitbx/vec3.h>
#include <cctbx/coordinates.h>
#include <cctbx/miller.h>

namespace cctbx {

  // forward declaration
  namespace sgtbx
  {
    class rot_mx;
    class change_of_basis_op;
  }

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
  inline af::shared<double>
  d_star_sq_as_stol_sq(af::const_ref<double> const &d_star_sq)
  {
    af::shared<double> result(
      d_star_sq.size(), af::init_functor_null<double>());
    for(std::size_t i=0;i<d_star_sq.size();i++) {
      result[i] = d_star_sq_as_stol_sq(d_star_sq[i]);
    }
    return result;
  }

  //! Conversion of d-spacing measures.
  inline double d_star_sq_as_two_stol(double d_star_sq)
  {
    return std::sqrt(d_star_sq);
  }

  //! Conversion of d-spacing measures.
  inline af::shared<double>
  d_star_sq_as_two_stol(af::const_ref<double> const &d_star_sq)
  {
    af::shared<double> result(
      d_star_sq.size(), af::init_functor_null<double>());
    for(std::size_t i=0;i<d_star_sq.size();i++) {
      result[i] = d_star_sq_as_two_stol(d_star_sq[i]);
    }
    return result;
  }

  //! Conversion of d-spacing measures.
  inline double d_star_sq_as_stol(double d_star_sq)
  {
    return std::sqrt(d_star_sq) * .5;
  }

  //! Conversion of d-spacing measures.
  inline af::shared<double>
  d_star_sq_as_stol(af::const_ref<double> const &d_star_sq)
  {
    af::shared<double> result(
      d_star_sq.size(), af::init_functor_null<double>());
    for(std::size_t i=0;i<d_star_sq.size();i++) {
      result[i] = d_star_sq_as_stol(d_star_sq[i]);
    }
    return result;
  }

  //! Conversion of d-spacing measures.
  inline double d_star_sq_as_d(double d_star_sq)
  {
    if (d_star_sq == 0.) return -1.;
    return 1. / std::sqrt(d_star_sq);
  }

  //! Conversion of d-spacing measures.
  inline af::shared<double>
  d_star_sq_as_d(af::const_ref<double> const &d_star_sq)
  {
    af::shared<double> result(
      d_star_sq.size(), af::init_functor_null<double>());
    for(std::size_t i=0;i<d_star_sq.size();i++) {
      result[i] = d_star_sq_as_d(d_star_sq[i]);
    }
    return result;
  }

  //! Conversion of d-spacing measures.
  inline double d_star_sq_as_two_theta(
    double d_star_sq, double wavelength, bool deg=false)
  {
    double sin_theta = d_star_sq_as_stol(d_star_sq) * wavelength;
    CCTBX_ASSERT(sin_theta <= 1.0);
    double result = 2. * std::asin(sin_theta);
    if (deg) return scitbx::rad_as_deg(result);
    return result;
  }

  //! Conversion of d-spacing measures.
  inline af::shared<double> d_star_sq_as_two_theta(
    af::const_ref<double> const &d_star_sq, double wavelength, bool deg=false)
  {
    af::shared<double> result(
      d_star_sq.size(), af::init_functor_null<double>());
    for(std::size_t i=0;i<d_star_sq.size();i++) {
      result[i] = d_star_sq_as_two_theta(d_star_sq[i], wavelength, deg);
    }
    return result;
  }

  //! Conversion of d-spacing measures.
  inline double stol_sq_as_d_star_sq(double stol_sq)
  {
    return stol_sq * 4;
  }

  //! Conversion of d-spacing measures.
  inline af::shared<double>
  stol_sq_as_d_star_sq(af::const_ref<double> const &stol_sq)
  {
    af::shared<double> result(
      stol_sq.size(), af::init_functor_null<double>());
    for(std::size_t i=0;i<stol_sq.size();i++) {
      result[i] = stol_sq_as_d_star_sq(stol_sq[i]);
    }
    return result;
  }

  //! Conversion of d-spacing measures.
  inline double two_stol_as_d_star_sq(double two_stol)
  {
    return std::pow(two_stol, 2);
  }

  //! Conversion of d-spacing measures.
  inline af::shared<double>
  two_stol_as_d_star_sq(af::const_ref<double> const &two_stol)
  {
    af::shared<double> result(
      two_stol.size(), af::init_functor_null<double>());
    for(std::size_t i=0;i<two_stol.size();i++) {
      result[i] = two_stol_as_d_star_sq(two_stol[i]);
    }
    return result;
  }

  //! Conversion of d-spacing measures.
  inline double stol_as_d_star_sq(double stol)
  {
    return std::pow(stol * 2, 2);
  }

  //! Conversion of d-spacing measures.
  inline af::shared<double>
  stol_as_d_star_sq(af::const_ref<double> const &stol)
  {
    af::shared<double> result(
      stol.size(), af::init_functor_null<double>());
    for(std::size_t i=0;i<stol.size();i++) {
      result[i] = stol_as_d_star_sq(stol[i]);
    }
    return result;
  }

  //! Conversion of d-spacing measures.
  inline double d_as_d_star_sq(double d)
  {
    if (d == 0.) return -1.;
    return 1. / std::pow(d, 2);
  }

  //! Conversion of d-spacing measures.
  inline af::shared<double>
  d_as_d_star_sq(af::const_ref<double> const &d)
  {
    af::shared<double> result(
      d.size(), af::init_functor_null<double>());
    for(std::size_t i=0;i<d.size();i++) {
      result[i] = d_as_d_star_sq(d[i]);
    }
    return result;
  }

  //! Conversion of d-spacing measures.
  inline double two_theta_as_d_star_sq(
    double two_theta, double wavelength, bool deg=false)
  {
    double theta = .5 * two_theta;
    if (deg) theta = scitbx::deg_as_rad(theta);
    double sin_theta = std::sin(theta);
    return stol_as_d_star_sq(sin_theta/wavelength);
  }

  //! Conversion of d-spacing measures.
  inline af::shared<double> two_theta_as_d_star_sq(
    af::const_ref<double> const &two_theta, double wavelength, bool deg=false)
  {
    af::shared<double> result(
      two_theta.size(), af::init_functor_null<double>());
    for(std::size_t i=0;i<two_theta.size();i++) {
      result[i] = two_theta_as_d_star_sq(two_theta[i], wavelength, deg);
    }
    return result;
  }

  //! Conversion of d-spacing measures.
  inline double two_theta_as_d(
    double two_theta, double wavelength, bool deg=false)
  {
    return d_star_sq_as_d(
      two_theta_as_d_star_sq(two_theta, wavelength, deg));
  }

  //! Conversion of d-spacing measures.
  inline af::shared<double> two_theta_as_d(
    af::const_ref<double> const &two_theta, double wavelength, bool deg=false)
  {
    af::shared<double> result(
      two_theta.size(), af::init_functor_null<double>());
    for(std::size_t i=0;i<two_theta.size();i++) {
      result[i] = two_theta_as_d(two_theta[i], wavelength, deg);
    }
    return result;
  }

  template <typename FloatType>
  FloatType
  mean_square_difference(
    scitbx::mat3<FloatType> const& a,
    scitbx::mat3<FloatType> const& b)
  {
    FloatType result = 0;
    for(std::size_t i=0;i<9;i++) {
      result += scitbx::fn::pow2(a[i] - b[i]);
    }
    return result;
  }

  //! Foadi, J., Evans, G. (2011). Acta Cryst. A67, 93-95.
  bool
  unit_cell_angles_are_feasible(
    scitbx::vec3<double> const& values_deg,
    double tolerance=1e-6);

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
          See also: paramters()
       */
      explicit
      unit_cell(af::small<double, 6> const& parameters);

      //! Constructor using parameters (a, b, c, alpha, beta, gamma).
      /*! See also: paramters()
       */
      explicit
      unit_cell(af::double6 const& parameters);

      //! Constructor using parameters derived from a metrical matrix.
      /*! The metrical matrix is defined as:
         <pre>
         ( a*a,            a*b*cos(gamma), a*c*cos(beta)  )
         ( a*b*cos(gamma), b*b,            b*c*cos(alpha) )
         ( a*c*cos(beta),  b*c*cos(alpha), c*c            )
         </pre>
          See also: metrical_matrix()
       */
      explicit
      unit_cell(uc_sym_mat3 const& metrical_matrix);

      //! Constructor using an orthogonalization matrix.
      /*! See also: orthogonalization_matrix()
       */
      explicit
      unit_cell(uc_mat3 const& orthogonalization_matrix);

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

      af::double6 const&
      d_volume_d_params() const { return d_volume_d_params_; }

      //! Corresponding reciprocal cell.
      unit_cell
      reciprocal() const;

      //! Length^2 of the longest lattice vector in the unit cell.
      double
      longest_vector_sq() const;

      //! Length^2 of the shortest lattice translation vector.
      /*! The fast_minimum_reduction is used to determine the shortest
          vector.
       */
      double
      shortest_vector_sq() const;

      //! Simple test for degenerated unit cell parameters.
      bool
      is_degenerate(double min_min_length_over_max_length=1.e-10,
                    double min_volume_over_min_length=1.e-5);

      //! Comparison of unit cell parameters.
      bool
      is_similar_to(unit_cell const& other,
                    double relative_length_tolerance=0.01,
                    double absolute_angle_tolerance=1.) const;

      //! Array of c_inv_r matrices compatible with change_basis.
      af::shared<scitbx::mat3<int> >
      similarity_transformations(
        unit_cell const& other,
        double relative_length_tolerance=0.02,
        double absolute_angle_tolerance=2.,
        int unimodular_generator_range=1) const;

      //! Matrix for the conversion of cartesian to fractional coordinates.
      /*! x(fractional) = matrix * x(cartesian). */
      uc_mat3 const& fractionalization_matrix() const { return frac_; }

      //! Matrix for the conversion of fractional to cartesian coordinates.
      /*! x(cartesian) = matrix * x(fractional). */
      uc_mat3 const& orthogonalization_matrix() const { return orth_; }

      //! Matrix for the conversion of grid indices to cartesian coordinates.
      /*! x(cartesian) = matrix * grid_index. */
      template <class IntType>
      uc_mat3
      grid_index_as_site_cart_matrix(
        scitbx::vec3<IntType> const& gridding) const
      {
        uc_mat3 result = orth_;
        for(unsigned i=0;i<3;i++) {
          CCTBX_ASSERT(gridding[i] > 0);
          double f = 1. / boost::numeric_cast<double>(gridding[i]);
          for(unsigned j=0;j<9;j+=3) result[i+j] *= f;
        }
        return result;
      }

      //! Conversion of cartesian coordinates to fractional coordinates.
      template <class FloatType>
      scitbx::vec3<FloatType>
      fractionalize(scitbx::vec3<FloatType> const& site_cart) const
      {
        // take advantage of the fact that frac_ is upper-triangular.
        return fractional<FloatType>(
            frac_[0] * site_cart[0]
          + frac_[1] * site_cart[1]
          + frac_[2] * site_cart[2],
            frac_[4] * site_cart[1]
          + frac_[5] * site_cart[2],
            frac_[8] * site_cart[2]);
      }

      //! Conversion of fractional coordinates to cartesian coordinates.
      template <class FloatType>
      scitbx::vec3<FloatType>
      orthogonalize(scitbx::vec3<FloatType> const& site_frac) const
      {
        // take advantage of the fact that orth_ is upper-triangular.
        return cartesian<FloatType>(
            orth_[0] * site_frac[0]
          + orth_[1] * site_frac[1]
          + orth_[2] * site_frac[2],
            orth_[4] * site_frac[1]
          + orth_[5] * site_frac[2],
            orth_[8] * site_frac[2]);
      }

      //! Conversion of cartesian coordinates to fractional coordinates.
      template <typename FloatType>
      af::shared<scitbx::vec3<FloatType> >
      fractionalize(
        af::const_ref<scitbx::vec3<FloatType> > const& sites_cart) const
      {
        af::shared<scitbx::vec3<FloatType> > result(
          sites_cart.size(),
          af::init_functor_null<scitbx::vec3<FloatType> >());
        const scitbx::vec3<FloatType>* si = sites_cart.begin();
        scitbx::vec3<FloatType>* ri = result.begin();
        for(std::size_t i=0;i<sites_cart.size();i++,ri++,si++) {
          // take advantage of the fact that frac_ is upper-triangular.
          (*ri)[0] = frac_[0] * (*si)[0]
                   + frac_[1] * (*si)[1]
                   + frac_[2] * (*si)[2];
          (*ri)[1] = frac_[4] * (*si)[1]
                   + frac_[5] * (*si)[2];
          (*ri)[2] = frac_[8] * (*si)[2];
        }
        return result;
      }

      //! Conversion of fractional coordinates to cartesian coordinates.
      template <typename FloatType>
      af::shared<scitbx::vec3<FloatType> >
      orthogonalize(
        af::const_ref<scitbx::vec3<FloatType> > const& sites_frac) const
      {
        af::shared<scitbx::vec3<FloatType> > result(
          sites_frac.size(),
          af::init_functor_null<scitbx::vec3<FloatType> >());
        const scitbx::vec3<FloatType>* si = sites_frac.begin();
        scitbx::vec3<FloatType>* ri = result.begin();
        for(std::size_t i=0;i<sites_frac.size();i++,ri++,si++) {
          // take advantage of the fact that orth_ is upper-triangular.
          (*ri)[0] = orth_[0] * (*si)[0]
                   + orth_[1] * (*si)[1]
                   + orth_[2] * (*si)[2];
          (*ri)[1] = orth_[4] * (*si)[1]
                   + orth_[5] * (*si)[2];
          (*ri)[2] = orth_[8] * (*si)[2];
        }
        return result;
      }

      //! Optimized fractionalization_matrix().transpose() * v
      /*! Not available in Python.
       */
      template <class FloatType>
      scitbx::vec3<FloatType>
      v_times_fractionalization_matrix_transpose(
        scitbx::vec3<FloatType> const& v) const
      {
        // take advantage of the fact that frac_ is upper-triangular.
        return fractional<FloatType>(
            frac_[0] * v[0],
            frac_[1] * v[0]
          + frac_[4] * v[1],
            frac_[2] * v[0]
          + frac_[5] * v[1]
          + frac_[8] * v[2]);
      }

      // ! gradient wrt Cartesian coordinates from gradient wrt fractional ones
      /*! The formula is \f$\nabla_{x_c} = F^T \nabla_{x_f}\f$ where \f$F\f$ is
        the fractionalisation matrix
        This is syntactic sugar for the long winded and cryptic
          v_times_fractionalization_matrix_transpose.
        */
      template <class FloatType>
      scitbx::vec3<FloatType>
      orthogonalize_gradient (scitbx::vec3<FloatType> const& v) const
      {
        return v_times_fractionalization_matrix_transpose(v);
      }

      // ! gradient wrt fractional coordinates from gradient wrt Cartesian ones
      /*! The formula is \f$\nabla_{x_f} = O^T \nabla_{x_c}\f$ where \f$O\f$ is
          the orthogonalisation matrix
      */
      template <class FloatType>
      scitbx::vec3<FloatType>
      fractionalize_gradient (const scitbx::vec3<FloatType>& g) const
      {
        const uc_mat3& O = orth_;
        return fractional<FloatType>(O[0]*g[0],
                                     O[1]*g[0] + O[4]*g[1],
                                     O[2]*g[0] + O[5]*g[1] + O[8]*g[2]);
      }

      /// The linear form transforming ADP u* into u_iso
      /** If
            U = { {u_0, u_3, u_4}, {u_3, u_1, u_5}, {u_4, u_5, u_2} }

          is the layout of u* as per scitbx::sym_mat3 and

            O ={ {o_0, o_1, o_2}, {0, o_4, o_5}, {0, 0, o_8} }

          is the layout of the orthogonalisation matrix as per scitbx::mat3,

            u_iso = 1/3 Tr( O U O^T )

                = o_0^2 u_0 + (o_1^2 + o_4^2) u_1 + (o_2^2 + o_5^2 + o_8^2) u_2
                  + 2 o_0 o_1 u_3 + 2 o_0 o_2 u_4 + 2 (o_1 o_2 + o_4 o_5) u_5
       */
      af::double6 const &u_star_to_u_iso_linear_form() const {
        return u_star_to_u_iso_linear_form_;
      }

      /// The linear map transforming ADP u* into u_cart
      /** If
            U = { {u_0, u_3, u_4}, {u_3, u_1, u_5}, {u_4, u_5, u_2} }

          is the layout of u* as per scitbx::sym_mat3 and

            O ={ {o_0, o_1, o_2}, {0, o_4, o_5}, {0, 0, o_8} }


          is the layout of the orthogonalisation matrix as per scitbx::mat3,

            u_cart = ( O u* O^T )

            can be alternatively be written as

            u_cart = L u*, where

            L ={ {o_0^2, o_1^2, o_2^2, 2 o_0 o_1, 2 o_0 o_2, 2 o_1 o_2},
                 {0, o_4^2, o_5^2, 0, 0, 2 o_4 o_5,
                 {0, 0, o_8^2, 0, 0, 0},
                 {0, o_1 o_4, o_2 o_5, o_0 o_4, o_0 o_5, o_2 o_4+o_1 o_5},
                 {0, 0, o_2 o_8, 0, o_0 o_8, o_1 o_8},
                 {0, 0, o_5 o_8, 0, 0, o_4 o_8}
               }
       */
      af::const_ref<double, af::mat_grid>
      u_star_to_u_cart_linear_map() const {
        return
        af::tiny_mat_const_ref<double, 6, 6>(u_star_to_u_cart_linear_map_);
      }

      /// The linear map transforming ADP u* into u_cif
      af::double6 const &u_star_to_u_cif_linear_map() const {
        return u_star_to_u_cif_linear_map_;
      }

      //! The gradient of the elements of the metrical matrix wrt the unit cell params
      af::const_ref<double, af::mat_grid>
      d_metrical_matrix_d_params() const {
        return
        af::tiny_mat_const_ref<double, 6, 6>(d_metrical_matrix_d_params_);
      }
      //! Length^2 of a vector of fractional coordinates.
      /*! Not available in Python.
       */
      template <class FloatType>
      FloatType
      length_sq(fractional<FloatType> const& site_frac) const
      {
        return orthogonalize(site_frac).length_sq();
      }

      //! Length of a vector of fractional coordinates.
      template <class FloatType>
      FloatType
      length(fractional<FloatType> const& site_frac) const
      {
        return std::sqrt(length_sq(site_frac));
      }

      //! Distance^2 between two vectors of fractional coordinates.
      /*! Not available in Python.
       */
      template <class FloatType>
      FloatType
      distance_sq(fractional<FloatType> const& site_frac_1,
                  fractional<FloatType> const& site_frac_2) const
      {
        return length_sq(fractional<FloatType>(site_frac_1 - site_frac_2));
      }

      //! Distance between two vectors of fractional coordinates.
      template <class FloatType>
      FloatType
      distance(fractional<FloatType> const& site_frac_1,
               fractional<FloatType> const& site_frac_2) const
      {
        return length(fractional<FloatType>(site_frac_1 - site_frac_2));
      }

      //! Angle in degrees formed by sites 1, 2 and 3 for fractional coordinates.
      template <class FloatType>
      boost::optional<FloatType>
      angle(fractional<FloatType> const& site_frac_1,
            fractional<FloatType> const& site_frac_2,
            fractional<FloatType> const& site_frac_3) const
      {
        cartesian<FloatType> vec_12 = orthogonalize(site_frac_1 - site_frac_2);
        cartesian<FloatType> vec_32 = orthogonalize(site_frac_3 - site_frac_2);
        FloatType length_12 = vec_12.length();
        FloatType length_32 = vec_32.length();
        if (length_12 == 0 || length_32 == 0) {
          return boost::optional<FloatType>();
        }
        const FloatType cos_angle = std::max(-1.,std::min(1.,
          (vec_12 * vec_32)/(length_12 * length_32)));
        return boost::optional<FloatType>(
          std::acos(cos_angle) / scitbx::constants::pi_180);
      }

      //! Dihedral angle in degrees formed by sites 1, 2, 3 and 4 for fractional coordinates.
      /*! angle = acos[(u x v).(v x w)/(|u x v||v x w|)]
          The sign of the angle is the sign of (u x v).w
       */
      template <class FloatType>
      boost::optional<FloatType>
      dihedral(fractional<FloatType> const& site_frac_1,
               fractional<FloatType> const& site_frac_2,
               fractional<FloatType> const& site_frac_3,
               fractional<FloatType> const& site_frac_4) const
      {
        cartesian<FloatType> u = orthogonalize(site_frac_1 - site_frac_2);
        cartesian<FloatType> v = orthogonalize(site_frac_3 - site_frac_2);
        cartesian<FloatType> w = orthogonalize(site_frac_3 - site_frac_4);
        cartesian<FloatType> u_cross_v = u.cross(v);
        cartesian<FloatType> v_cross_w = v.cross(w);
        FloatType norm_u_cross_v = u_cross_v.length_sq();
        if (norm_u_cross_v == 0) return boost::optional<FloatType>();
        FloatType norm_v_cross_w = v_cross_w.length_sq();
        if (norm_v_cross_w == 0) return boost::optional<FloatType>();
        FloatType cos_angle = std::max(-1.,std::min(1.,
          u_cross_v * v_cross_w / std::sqrt(norm_u_cross_v * norm_v_cross_w)));
        FloatType angle = std::acos(cos_angle) / scitbx::constants::pi_180;
        if (u_cross_v * w < 0) {
          angle *= -1;
        }
        return boost::optional<FloatType>(angle);
      }

      /*! \brief Shortest length^2 of a vector of fractional coordinates
          under application of periodicity. */
      /*! Not available in Python.
       */
      template <class FloatType>
      FloatType
      mod_short_length_sq(fractional<FloatType> const& site_frac) const
      {
        return length_sq(site_frac.mod_short());
      }

      /*! \brief Shortest length of a vector of fractional coordinates
          under application of periodicity. */
      template <class FloatType>
      FloatType
      mod_short_length(fractional<FloatType> const& site_frac) const
      {
        return std::sqrt(mod_short_length_sq(site_frac));
      }

      /*! \brief Shortest distance^2 between two vectors of fractional
          coordinates under application of periodicity. */
      /*! Not available in Python.
       */
      template <class FloatType>
      FloatType
      mod_short_distance_sq(fractional<FloatType> const& site_frac_1,
                            fractional<FloatType> const& site_frac_2) const
      {
        return mod_short_length_sq(
          fractional<FloatType>(site_frac_1 - site_frac_2));
      }

      /*! \brief Shortest distance between two vectors of fractional
          coordinates under application of periodicity. */
      template <class FloatType>
      FloatType
      mod_short_distance(fractional<FloatType> const& site_frac_1,
                         fractional<FloatType> const& site_frac_2) const
      {
        return std::sqrt(mod_short_distance_sq(site_frac_1, site_frac_2));
      }

      /*! \brief Shortest distance^2 between all sites in site_frac_1
          and site_frac_2 under application of periodicity.
       */
      /*! Not available in Python.
       */
      template <class FloatType>
      FloatType
      min_mod_short_distance_sq(
        af::const_ref<scitbx::vec3<FloatType> > const& sites_frac_1,
        fractional<FloatType> const& site_frac_2) const
      {
        FloatType
          result = mod_short_distance_sq(
            fractional<FloatType>(sites_frac_1[0]), site_frac_2);
        for(std::size_t i=1;i<sites_frac_1.size();i++) {
          scitbx::math::update_min(
            result, mod_short_distance_sq(
              fractional<FloatType>(sites_frac_1[i]), site_frac_2));
        }
        return result;
      }

      /*! \brief Shortest distance between all sites in sites_frac_1
          and site_frac_2 under application of periodicity.
       */
      template <class FloatType>
      FloatType
      min_mod_short_distance(
        af::const_ref<scitbx::vec3<FloatType> > const& sites_frac_1,
        fractional<FloatType> const& site_frac_2) const
      {
        return std::sqrt(min_mod_short_distance_sq(sites_frac_1, site_frac_2));
      }

      //! Transformation (change-of-basis) of unit cell parameters.
      /*! c_inv_r is the inverse of the 3x3 change-of-basis matrix
          that transforms coordinates in the old basis system to
          coodinates in the new basis system.
       */
      unit_cell
      change_basis(uc_mat3 const& c_inv_r, double r_den=1.) const;

      //! Transformation (change-of-basis) of unit cell parameters.
      /*! r is the inverse of the 3x3 change-of-basis matrix
          that transforms coordinates in the old basis system to
          coodinates in the new basis system.
          <p>
          Not available in Python.
       */
      unit_cell
      change_basis(sgtbx::rot_mx const& c_inv_r) const;

      //! Transformation (change-of-basis) of unit cell parameters.
      unit_cell
      change_basis(sgtbx::change_of_basis_op const& cb_op) const;

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
      d_star_sq(miller::index<NumType> const& miller_index) const
      {
        return
            (miller_index[0] * miller_index[0]) * r_metr_mx_[0]
          + (miller_index[1] * miller_index[1]) * r_metr_mx_[1]
          + (miller_index[2] * miller_index[2]) * r_metr_mx_[2]
          + (2 * miller_index[0] * miller_index[1]) * r_metr_mx_[3]
          + (2 * miller_index[0] * miller_index[2]) * r_metr_mx_[4]
          + (2 * miller_index[1] * miller_index[2]) * r_metr_mx_[5];
      }

      //! d-spacing measure d_star_sq = 1/d^2 = (2*sin(theta)/lambda)^2.
      template <typename NumType>
      af::shared<double>
      d_star_sq(
        af::const_ref<miller::index<NumType> > const& miller_indices) const
      {
        af::shared<double> result(
          miller_indices.size(), af::init_functor_null<double>());
        for(std::size_t i=0;i<miller_indices.size();i++) {
          result[i] = d_star_sq(miller_indices[i]);
        }
        return result;
      }

      //! Maximum d_star_sq for given list of Miller indices.
      template <typename NumType>
      double
      max_d_star_sq(
        af::const_ref<miller::index<NumType> > const& miller_indices) const
      {
        double result = 0;
        for(std::size_t i=0;i<miller_indices.size();i++) {
          scitbx::math::update_max(result, d_star_sq(miller_indices[i]));
        }
        return result;
      }

      //! Minimum and maximum d_star_sq for given list of Miller indices.
      template <typename NumType>
      af::double2
      min_max_d_star_sq(
        af::const_ref<miller::index<NumType> > const& miller_indices) const
      {
        af::double2 result(0, 0);
        if (miller_indices.size()) {
          result.fill(d_star_sq(miller_indices[0]));
          for(std::size_t i=1;i<miller_indices.size();i++) {
            double q = d_star_sq(miller_indices[i]);
            scitbx::math::update_min(result[0], q);
            scitbx::math::update_max(result[1], q);
          }
        }
        return result;
      }

      //! d-spacing measure (sin(theta)/lambda)^2 = d_star_sq/4.
      template <typename NumType>
      double
      stol_sq(miller::index<NumType> const& miller_index) const
      {
        return d_star_sq_as_stol_sq(d_star_sq(miller_index));
      }

      //! d-spacing measure (sin(theta)/lambda)^2 = d_star_sq/4.
      template <typename NumType>
      af::shared<double>
      stol_sq(
        af::const_ref<miller::index<NumType> > const& miller_indices) const
      {
        af::shared<double> result(
          miller_indices.size(), af::init_functor_null<double>());
        for(std::size_t i=0;i<miller_indices.size();i++) {
          result[i] = stol_sq(miller_indices[i]);
        }
        return result;
      }

      //! d-spacing measure 2*sin(theta)/lambda = 1/d = sqrt(d_star_sq).
      template <typename NumType>
      double
      two_stol(miller::index<NumType> const& miller_index) const
      {
        return d_star_sq_as_two_stol(d_star_sq(miller_index));
      }

      //! d-spacing measure 2*sin(theta)/lambda = 1/d = sqrt(d_star_sq).
      template <typename NumType>
      af::shared<double>
      two_stol(
        af::const_ref<miller::index<NumType> > const& miller_indices) const
      {
        af::shared<double> result(
          miller_indices.size(), af::init_functor_null<double>());
        for(std::size_t i=0;i<miller_indices.size();i++) {
          result[i] = two_stol(miller_indices[i]);
        }
        return result;
      }

      //! d-spacing measure sin(theta)/lambda = 1/(2*d) = sqrt(d_star_sq)/2.
      template <typename NumType>
      double
      stol(miller::index<NumType> const& miller_index) const
      {
        return d_star_sq_as_stol(d_star_sq(miller_index));
      }

      //! d-spacing measure sin(theta)/lambda = 1/(2*d) = sqrt(d_star_sq)/2.
      template <typename NumType>
      af::shared<double>
      stol(af::const_ref<miller::index<NumType> > const& miller_indices) const
      {
        af::shared<double> result(
          miller_indices.size(), af::init_functor_null<double>());
        for(std::size_t i=0;i<miller_indices.size();i++) {
          result[i] = stol(miller_indices[i]);
        }
        return result;
      }

      //! d-spacing measure d = 1/(2*sin(theta)/lambda).
      template <typename NumType>
      double
      d(miller::index<NumType> const& miller_index) const
      {
        return d_star_sq_as_d(d_star_sq(miller_index));
      }

      //! d-spacing measure d = 1/(2*sin(theta)/lambda).
      template <typename NumType>
      double
      d_frac(scitbx::vec3<NumType> const& miller_index) const
      {
        return d_star_sq_as_d(d_star_sq(miller::index<NumType>(miller_index)));
      }

      //! d-spacing measure d = 1/(2*sin(theta)/lambda).
      template <typename NumType>
      af::shared<double>
      d(af::const_ref<miller::index<NumType> > const& miller_indices) const
      {
        af::shared<double> result(
          miller_indices.size(), af::init_functor_null<double>());
        for(std::size_t i=0;i<miller_indices.size();i++) {
          result[i] = d(miller_indices[i]);
        }
        return result;
      }

      //! Diffraction angle 2-theta, given wavelength.
      template <typename NumType>
      double
      two_theta(
        miller::index<NumType> const& miller_index,
        double wavelength,
        bool deg=false) const
      {
        return d_star_sq_as_two_theta(
          d_star_sq(miller_index), wavelength, deg);
      }

      //! Diffraction angle 2-theta, given wavelength.
      template <typename NumType>
      af::shared<double>
      two_theta(
        af::const_ref<miller::index<NumType> > const& miller_indices,
        double wavelength,
        bool deg=false) const
      {
        af::shared<double> result(
          miller_indices.size(), af::init_functor_null<double>());
        for(std::size_t i=0;i<miller_indices.size();i++) {
          result[i] = two_theta(miller_indices[i], wavelength, deg);
        }
        return result;
      }

      //! Squared sine of the diffraction angle 2-theta, given wavelength.
      /*! The result is always positive, using  trigonometric identities:
          sin(2a) = 2*sin(a)*cos(a) and cos^2(a) = 1-sin^2(a)
      */
      template <typename NumType>
      double
      sin_sq_two_theta(
        miller::index<NumType> const& miller_index,
        double wavelength) const
      {
        const double x = d_star_sq_as_stol_sq(d_star_sq(miller_index)) *
          std::pow(wavelength, 2.0);
        return std::max(0.0, 4*x*(1-x));
      }

      //! Squared sine of the diffraction angle 2-theta, given wavelength.
      template <typename NumType>
      af::shared<double>
      sin_sq_two_theta(
        af::const_ref<miller::index<NumType> > const& miller_indices,
        double wavelength) const
      {
        af::shared<double> result(
          miller_indices.size(), af::init_functor_null<double>());
        for(std::size_t i=0;i<miller_indices.size();i++) {
          result[i] = sin_sq_two_theta(miller_indices[i], wavelength);
        }
        return result;
      }

      //! Sine of the diffraction angle 2-theta, given wavelength.
      template <typename NumType>
      double
      sin_two_theta(
        miller::index<NumType> const& miller_index,
        double wavelength) const
      {
        const double x = d_star_sq_as_stol_sq(d_star_sq(miller_index)) *
          std::pow(wavelength, 2.0);
        return 2*std::sqrt(std::max(0.0, x*(1-x)));
      }

      //! Sine of the diffraction angle 2-theta, given wavelength.
      template <typename NumType>
      af::shared<double>
      sin_two_theta(
        af::const_ref<miller::index<NumType> > const& miller_indices,
        double wavelength) const
      {
        af::shared<double> result(
          miller_indices.size(), af::init_functor_null<double>());
        for(std::size_t i=0;i<miller_indices.size();i++) {
          result[i] = sin_two_theta(miller_indices[i], wavelength);
        }
        return result;
      }

      template <typename NumType>
      scitbx::vec3<double>
      reciprocal_space_vector(
        miller::index<NumType> const& miller_index) const
      {
        uc_mat3 const& frac_mat = fractionalization_matrix();
        scitbx::vec3<double> rcv = miller_index * frac_mat;
        return rcv;
      }

      template <typename NumType>
      af::shared<scitbx::vec3<double> >
      reciprocal_space_vector(
        af::const_ref<miller::index<NumType> > const& miller_indices) const
      {
        af::shared<scitbx::vec3<double> > result(
          miller_indices.size(),
          af::init_functor_null<scitbx::vec3<double> >());
        for (std::size_t i = 0; i < miller_indices.size(); i++) {
          result[i] = reciprocal_space_vector(miller_indices[i]);
        }
        return result;
      }

      //! Simple measure for the similarity of two unit cells.
      /*! The result is the mean of the squared differences between
          basis vectors. The basis vectors are defined by the columns
          of this->orthogonalization_matrix() and
          other.orthogonalization_matrix().
       */
      double
      bases_mean_square_difference(unit_cell const& other) const
      {
        return mean_square_difference(orth_, other.orth_) / 3;
      }

      /*! \brief Sort order for alternative settings of otherwise
          identical orthorhombic cells.
       */
      /*! The cell with the smaller a-parameter is preferred.
          If the a-parameters are equal, the cell with the
          smaller b-parameter is preferred.
          If the a-parameters and b-parameters are equal, the cell
          with the smaller c-parameter is preferred.

          The return value is -1 if this is the preferred cell,
          1 if other is the preferred cell, and 0 otherwise.

          The result is meaningless if the two unit cells are
          not alternative settings; i.e. the two cells must be
          related to each other by a transformation matrix
          with determinant 1.
       */
      int
      compare_orthorhombic(
        const unit_cell& other) const;

      /*! \brief Sort order for alternative settings of otherwise
          identical monoclinic cells.
       */
      /*! unique_axis must be one of:
            0 for a-unique,
            1 for b-unique,
            2 for c-unique.

          If the absolute value of the difference between the
          monoclinic angels is less than angular_tolerance
          compare_orthorhombic() is used to determine the return value
          (-1, 1, or 0). Otherwise the return value is determined
          by evaluating the monoclinic angels. For exact details
          refer to the source code of the implementation
          ($CCTBX_DIST/cctbx/uctbx/uctbx.cpp).

          The result is meaningless if the two unit cells are
          not alternative settings; i.e. the two cells must be
          related to each other by a transformation matrix
          with determinant 1.
       */
      int
      compare_monoclinic(
        const unit_cell& other,
        unsigned unique_axis,
        double angular_tolerance) const;

      sgtbx::change_of_basis_op const&
      change_of_basis_op_for_best_monoclinic_beta() const;

    protected:
      void init_volume();
      void init_reciprocal();
      void init_metrical_matrices();
      void init_orth_and_frac_matrices();
      void init_tensor_rank_2_orth_and_frac_linear_maps();
      void initialize();

      af::double6 params_;
      af::double3 sin_ang_;
      af::double3 cos_ang_;
      double volume_;
      af::double6 d_volume_d_params_;
      uc_sym_mat3 metr_mx_;
      af::double6 r_params_;
      af::double3 r_sin_ang_;
      af::double3 r_cos_ang_;
      uc_sym_mat3 r_metr_mx_;
      uc_mat3 frac_;
      uc_mat3 orth_;
      af::double6 u_star_to_u_iso_linear_form_;
      af::double6 u_star_to_u_cif_linear_map_;
      double u_star_to_u_cart_linear_map_[6*6];
      double d_metrical_matrix_d_params_[6*6];

      mutable double longest_vector_sq_;
      mutable double shortest_vector_sq_;

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

  struct distance_mod_1
  {
    fractional<> diff_raw;
    fractional<> diff_mod;
    double dist_sq;

    distance_mod_1() {}

    distance_mod_1(
      uctbx::unit_cell const& unit_cell,
      fractional<> const& site_frac_1,
      fractional<> const& site_frac_2)
    :
      diff_raw(site_frac_1 - site_frac_2),
      diff_mod(diff_raw.mod_short()),
      dist_sq(unit_cell.length_sq(diff_mod))
    {}

    scitbx::vec3<int>
    unit_shifts() const
    {
      return fractional<>(diff_mod - diff_raw).unit_shifts();
    }
  };

}} // namespace cctbx::uctbx

#endif // CCTBX_UCTBX_H
