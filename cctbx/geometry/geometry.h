#ifndef CCTBX_GEOMETRY_H
#define CCTBX_GEOMETRY_H

#include <cctbx/import_scitbx_af.h>
#include <cctbx/uctbx.h>
#include <cctbx/sgtbx/rt_mx.h>
#include <scitbx/array_family/accessors/packed_matrix.h>
#include <scitbx/matrix/matrix_vector_operations.h>
#include <tbxx/optional_copy.hpp>

#include <vector>

namespace cctbx { namespace geometry {

  namespace af=scitbx::af;
  using tbxx::optional_container;

  namespace detail {

    template<typename af_tiny_t, typename FloatType>
    FloatType
    variance_impl(
      af_tiny_t const &grads,
      af::const_ref<FloatType, af::packed_u_accessor> const &covariance_matrix)
    {
      FloatType result = 0;
      for (std::size_t i_seq=0; i_seq<grads.size(); i_seq++) {
        for (std::size_t j_seq=i_seq; j_seq<grads.size(); j_seq++) {
          for (std::size_t i=0; i<3; i++) {
            for (std::size_t j=0; j<3; j++) {
              if (i_seq==j_seq && j<i) { continue; }
              FloatType tmp = grads[i_seq][i] * grads[j_seq][j] *
                              covariance_matrix(i_seq*3+i, j_seq*3+j);
              if ( !(i_seq==j_seq && i==j) ) tmp *= 2.;
              result += tmp;
            }
          }
        }
      }
      return result;
    }
  } // namespace detail

  template<typename FloatType>
  class distance
  {
  public:
    distance() {}

    distance(af::tiny<scitbx::vec3<FloatType>, 2> const &sites_)
      :
    sites(sites_)
    {
      init_distance_model();
    }

    scitbx::vec3<FloatType>
    d_distance_d_site_0(FloatType epsilon=1.e-100) const
    {
      if (distance_model < epsilon) return scitbx::vec3<FloatType>(0,0,0);
      return (sites[1] - sites[0]) / distance_model;
    }

    af::tiny<scitbx::vec3<FloatType>, 2>
    d_distance_d_sites(FloatType epsilon=1.e-100) const
    {
      af::tiny<scitbx::vec3<FloatType>, 2> result;
      result[0] = d_distance_d_site_0();
      result[1] = -result[0];
      return result;
    }

    // The gradient of the distance wrt the elements of the metrical matrix
    scitbx::sym_mat3<FloatType>
    d_distance_d_metrical_matrix(cctbx::uctbx::unit_cell const &unit_cell) const
    {
      scitbx::vec3<FloatType> vec_frac = unit_cell.fractionalize(sites[0]-sites[1]);
      scitbx::sym_mat3<FloatType> result;
      FloatType one_over_distance = 1./distance_model;
      result[0] = vec_frac[0] * vec_frac[0] * one_over_distance * 0.5;
      result[1] = vec_frac[1] * vec_frac[1] * one_over_distance * 0.5;
      result[2] = vec_frac[2] * vec_frac[2] * one_over_distance * 0.5;
      result[3] = vec_frac[0] * vec_frac[1] * one_over_distance;
      result[4] = vec_frac[0] * vec_frac[2] * one_over_distance;
      result[5] = vec_frac[1] * vec_frac[2] * one_over_distance;
      return result;
    }

    // The gradient of the distance wrt the elements of the unit cell parameters
    scitbx::sym_mat3<FloatType>
    d_distance_d_cell_params(cctbx::uctbx::unit_cell const &unit_cell) const
    {
      scitbx::sym_mat3<FloatType> result;
      scitbx::sym_mat3<FloatType>
        d_d_mm = d_distance_d_metrical_matrix(unit_cell);
      af::versa<FloatType, af::c_grid<2> >
        d_mm_d_params = unit_cell.d_metrical_matrix_d_params();
      scitbx::matrix::matrix_transposed_vector(
        6, 6, d_mm_d_params.begin(), d_d_mm.begin(), result.begin());
      return result;
    }

    // The variance of the distance taking into account only the errors in the sites
    FloatType
    variance(
      af::const_ref<FloatType, af::packed_u_accessor> const &covariance_matrix,
      cctbx::uctbx::unit_cell const &unit_cell,
      sgtbx::rt_mx const &rt_mx_ji) const
    {
      CCTBX_ASSERT(covariance_matrix.size() == 21);
      af::tiny<scitbx::vec3<FloatType>, 2> grads;
      grads[0] = d_distance_d_site_0();
      grads[1] = -grads[0];
      if (!rt_mx_ji.is_unit_mx()) {
        scitbx::mat3<double> r_inv_cart
          =   unit_cell.orthogonalization_matrix()
            * rt_mx_ji.r().inverse().as_double()
            * unit_cell.fractionalization_matrix();
        grads[1] = r_inv_cart * grads[1];
      }
      return detail::variance_impl(grads, covariance_matrix);
    }

    /*! The variance of the distance taking into account errors in the sites
        and errors in the unit cell parameters.

        Under the assumption that the errors in the sites are uncorrelated with
        the errors in the unit cell parameters, then

          sigma^2(f) = sigma^2(f,sites) + sigma^2(f,cell)
     */
    FloatType
    variance(
      af::const_ref<FloatType, af::packed_u_accessor> const &covariance_matrix,
      af::const_ref<FloatType, af::packed_u_accessor> const &
        cell_covariance_matrix,
      cctbx::uctbx::unit_cell const &unit_cell,
      sgtbx::rt_mx const &rt_mx_ji) const
    {
      CCTBX_ASSERT(cell_covariance_matrix.size() == 21);
      FloatType var = variance(covariance_matrix, unit_cell, rt_mx_ji);
      scitbx::sym_mat3<FloatType>
        d_distance_d_cell = d_distance_d_cell_params(unit_cell);
      var += scitbx::matrix::quadratic_form_packed_u(
        6, cell_covariance_matrix.begin(), d_distance_d_cell.begin());
      return var;
    }

    //! Cartesian coordinates of bonded sites.
    af::tiny<scitbx::vec3<FloatType>, 2> sites;
    //! Distance between sites.
    FloatType distance_model;

  protected:
    void
      init_distance_model()
    {
      distance_model = (sites[0] - sites[1]).length();
    }

  };

  template<typename FloatType>
  class angle
  {
  public:
    angle() {}

    angle(af::tiny<scitbx::vec3<FloatType>, 3> const &sites_)
      :
    sites(sites_)
    {
      init_angle_model();
    }

    //! Gradient of the angle with respect to the three sites.
    /*! The formula for the gradients is singular at delta = 0
        and delta = 180. However, the gradients converge to zero
        near these singularities. To avoid numerical problems, the
        gradients and curvatures are set to zero exactly if the
        intermediate result sqrt(1-cos(angle_model)**2) < epsilon.

        See also:
          http://salilab.org/modeller/manual/manual.html,
          "Features and their derivatives"
     */
    af::tiny<scitbx::vec3<FloatType>, 3>
    d_angle_d_sites(FloatType epsilon=1.e-100) const
    {
      af::tiny<scitbx::vec3<FloatType>, 3> grads;
      if (!have_angle_model) {
        std::fill_n(grads.begin(), 3U, scitbx::vec3<FloatType>(0,0,0));
      }
      else {
        FloatType
        sin_angle_model = std::sqrt(1-scitbx::fn::pow2(cos_angle_model));
        if (sin_angle_model < epsilon) {
          std::fill_n(grads.begin(), 3U, scitbx::vec3<FloatType>(0,0,0));
        }
        else {
          using scitbx::constants::pi_180;
          scitbx::vec3<FloatType> d_angle_d_site0, d_angle_d_site2;
          d_angle_d_site0 = (d_01_unit * cos_angle_model - d_21_unit) /
                            (sin_angle_model * d_01_abs);
          grads[0] = -d_angle_d_site0 / pi_180;
          d_angle_d_site2 = (d_21_unit * cos_angle_model - d_01_unit) /
                            (sin_angle_model * d_21_abs);
          grads[2] = -d_angle_d_site2 / pi_180;
          grads[1] = -(grads[0] + grads[2]);
        }
      }
      return grads;
    }

    // The gradient of the angle wrt the elements of the metrical matrix
    scitbx::sym_mat3<FloatType>
    d_angle_d_metrical_matrix(cctbx::uctbx::unit_cell const &unit_cell) const
    {
      scitbx::vec3<FloatType> d_01_frac = unit_cell.fractionalize(d_01);
      scitbx::vec3<FloatType> d_21_frac = unit_cell.fractionalize(d_21);
      scitbx::sym_mat3<FloatType> result;
      FloatType overall_factor
        = 1./std::sqrt(1-scitbx::fn::pow2(cos_angle_model));
      FloatType factor0 = cos_angle_model/(d_01_abs*d_01_abs);
      FloatType factor1 = 1./(d_01_abs*d_21_abs);
      FloatType factor2 = cos_angle_model/(d_21_abs*d_21_abs);
      for (std::size_t i=0;i<3;i++) {
        result[i] = 0.5 * overall_factor * (
            scitbx::fn::pow2(d_01_frac[i])* factor0
          - d_01_frac[i] * d_21_frac[i] * 2 * factor1
          + scitbx::fn::pow2(d_21_frac[i]) * factor2);
      }
      result[3] = overall_factor * (
           d_01_frac[0] * d_01_frac[1] * factor0
        - (d_01_frac[0] * d_21_frac[1] + d_01_frac[1] * d_21_frac[0]) * factor1
        +  d_21_frac[0] * d_21_frac[1] * factor2);
      result[4] = overall_factor * (
           d_01_frac[0] * d_01_frac[2] * factor0
        - (d_01_frac[0] * d_21_frac[2] + d_01_frac[2] * d_21_frac[0]) * factor1
        +  d_21_frac[0] * d_21_frac[2] * factor2);
      result[5] = overall_factor * (
           d_01_frac[1] * d_01_frac[2] * factor0
        - (d_01_frac[1] * d_21_frac[2] + d_01_frac[2] * d_21_frac[1]) * factor1
        +  d_21_frac[1] * d_21_frac[2] * factor2);
      return result;
    }

    // The gradient of the angle wrt the elements of the unit cell parameters
    scitbx::sym_mat3<FloatType>
    d_angle_d_cell_params(cctbx::uctbx::unit_cell const &unit_cell) const
    {
      scitbx::sym_mat3<FloatType> result;
      scitbx::sym_mat3<FloatType> d_d_mm = d_angle_d_metrical_matrix(unit_cell);
      af::versa<FloatType, af::c_grid<2> >
        d_mm_d_params = unit_cell.d_metrical_matrix_d_params();
      scitbx::matrix::matrix_transposed_vector(
        6, 6, d_mm_d_params.begin(), d_d_mm.begin(), result.begin());
      return result;
    }

    // The variance of the angle taking into account only the errors in the sites
    FloatType
    variance(
      af::const_ref<FloatType, af::packed_u_accessor> const &covariance_matrix,
      cctbx::uctbx::unit_cell const &unit_cell,
      optional_container<af::shared<sgtbx::rt_mx> > const &sym_ops) const
    {
      CCTBX_ASSERT(covariance_matrix.size() == 45);
      af::tiny<scitbx::vec3<FloatType>, 3> grads = d_angle_d_sites();
      for (std::size_t i=0; i<3; i++) {
        if (sym_ops && !sym_ops[i].is_unit_mx()) {
          scitbx::mat3<double> r_inv_cart
            =   unit_cell.orthogonalization_matrix()
              * sym_ops[i].r().inverse().as_double()
              * unit_cell.fractionalization_matrix();
          grads[i] = r_inv_cart * grads[i];
        }
      }
      return detail::variance_impl(grads, covariance_matrix);
    }

    /*! The variance of the angle taking into account errors in the sites
        and errors in the unit cell parameters.

        Under the assumption that the errors in the sites are uncorrelated with
        the errors in the unit cell parameters, then

          sigma^2(f) = sigma^2(f,sites) + sigma^2(f,cell)
     */
    FloatType
    variance(
      af::const_ref<FloatType, af::packed_u_accessor> const &covariance_matrix,
      af::const_ref<FloatType, af::packed_u_accessor> const &
        cell_covariance_matrix,
      cctbx::uctbx::unit_cell const &unit_cell,
      optional_container<af::shared<sgtbx::rt_mx> > const &sym_ops) const
    {
      CCTBX_ASSERT(cell_covariance_matrix.size() == 21);
      FloatType var = variance(covariance_matrix, unit_cell, sym_ops);
      scitbx::sym_mat3<FloatType>
        d_angle_d_cell = d_angle_d_cell_params(unit_cell);
      var += scitbx::matrix::quadratic_form_packed_u(
        6, cell_covariance_matrix.begin(), d_angle_d_cell.begin());
      return var;
    }

    //! Cartesian coordinates of sites forming the angle.
    af::tiny<scitbx::vec3<FloatType>, 3> sites;
    //! false in singular situations.
    bool have_angle_model;
    //! Value of angle formed by the sites.
    FloatType angle_model;

  protected:
    FloatType d_01_abs;
    FloatType d_21_abs;
    scitbx::vec3<FloatType> d_01;
    scitbx::vec3<FloatType> d_21;
    scitbx::vec3<FloatType> d_01_unit;
    scitbx::vec3<FloatType> d_21_unit;
    FloatType cos_angle_model;

    void
    init_angle_model()
    {
      have_angle_model = false;
      d_01_abs = 0;
      d_21_abs = 0;
      d_01.fill(0);
      d_21.fill(0);
      d_01_unit.fill(0);
      d_21_unit.fill(0);
      cos_angle_model = -9;
      d_01 = sites[0] - sites[1];
      d_01_abs = d_01.length();
      if (d_01_abs > 0) {
        d_21 = sites[2] - sites[1];
        d_21_abs = d_21.length();
        if (d_21_abs > 0) {
          d_01_unit = d_01 / d_01_abs;
          d_21_unit = d_21 / d_21_abs;
          cos_angle_model = std::max(-1.,std::min(1.,d_01_unit*d_21_unit));
          angle_model = std::acos(cos_angle_model)
                      / scitbx::constants::pi_180;
          have_angle_model = true;
        }
      }
    }
  };

}} //cctbx::geometry

#endif
