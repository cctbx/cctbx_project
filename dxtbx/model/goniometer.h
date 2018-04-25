/*
 * goniometer.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DXTBX_MODEL_GONIOMETER_H
#define DXTBX_MODEL_GONIOMETER_H

#include <iostream>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/math/r3_rotation.h>
#include <scitbx/constants.h>
#include <scitbx/array_family/simple_io.h>
#include <scitbx/array_family/simple_tiny_io.h>
#include <scitbx/array_family/shared.h>
#include <dxtbx/error.h>
#include "model_helpers.h"

namespace dxtbx { namespace model {

  using scitbx::vec3;
  using scitbx::mat3;
  using scitbx::constants::pi;
  using scitbx::math::r3_rotation::axis_and_angle_as_matrix;

  /** A goniometer base class */
  class GoniometerBase {};

  /**
   * A class to represent the rotation axis for a standard rotation
   * geometry diffraction data set.
   *
   * The rotation axis assumed to have its origin at the origin of the
   * laboratory coordinate system. The rotation axis vector is normalized
   * when set using either the constructor or the rotation axis 'setter'
   *
   * The fixed rotation matrix, F, represents the additional fixed rotation of
   * the sample attached to the rotating axis - for example the effects of
   * kappa and phi for an omega scan on a kappa goniometer. A second rotation
   * matrix, S, may be applied which is used to align the rotation axis,
   * for example the omega & kappa axes for a phi scan. These rotation matrices
   * should be applied as:
   *
   *  A = [S][R][F][U][B]
   *
   * in a standard orientation matrix.
   *
   * In the document detailing the conventions used:
   *    rotation_axis -> m2
   *    fixed_rotation_matrix -> *unspecified*
   */
  class Goniometer : public GoniometerBase {

  public:

    /**
     * Initialise the goniometer. The fixed rotation matrix is set to the
     * identity matrix and the rotation axis is set to the x axis
     * @param rotation_axis The goniometer rotation axis
     */
    Goniometer()
      : rotation_axis_(1.0, 0.0, 0.0),
        fixed_rotation_(
          1.0, 0.0, 0.0,
          0.0, 1.0, 0.0,
          0.0, 0.0, 1.0),
        setting_rotation_(
          1.0, 0.0, 0.0,
          0.0, 1.0, 0.0,
          0.0, 0.0, 1.0) {}

    /**
     * Initialise the goniometer. The fixed rotation matrix is set to the
     * identity matrix.
     * @param rotation_axis The goniometer rotation axis
     */
    Goniometer(vec3 <double> rotation_axis)
      : fixed_rotation_(
          1.0, 0.0, 0.0,
          0.0, 1.0, 0.0,
          0.0, 0.0, 1.0),
        setting_rotation_(
          1.0, 0.0, 0.0,
          0.0, 1.0, 0.0,
          0.0, 0.0, 1.0) {
      DXTBX_ASSERT(rotation_axis.length() > 0);
      rotation_axis_ = rotation_axis.normalize();
    }

    /**
     * Initialise the goniometer.
     * @param rotation_axis The goniometer rotation axis
     * @param fixed_rotation The additional fixed rotation of the sample
     */
    Goniometer(vec3 <double> rotation_axis,
               mat3 <double> fixed_rotation)
      : fixed_rotation_(fixed_rotation),
        setting_rotation_(
          1.0, 0.0, 0.0,
          0.0, 1.0, 0.0,
          0.0, 0.0, 1.0) {
      DXTBX_ASSERT(rotation_axis.length() > 0);
      rotation_axis_ = rotation_axis.normalize();
    }

    /**
     * Initialise the goniometer.
     * @param rotation_axis The goniometer rotation axis at datum
     * @param fixed_rotation The additional fixed rotation of the sample
     * @param setting_rotation The additional setting rotation of the axis
     */
    Goniometer(vec3 <double> rotation_axis,
               mat3 <double> fixed_rotation,
               mat3 <double> setting_rotation)
      : fixed_rotation_(fixed_rotation),
        setting_rotation_(setting_rotation) {
      DXTBX_ASSERT(rotation_axis.length() > 0);
      rotation_axis_ = rotation_axis.normalize();
    }

    /** Virtual destructor */
    virtual ~Goniometer() {}

    /** Get the rotation axis */
    vec3 <double> get_rotation_axis() const {
      return setting_rotation_ * rotation_axis_;
    }

    /** Get the rotation axis at datum */
    vec3 <double> get_rotation_axis_datum() const {
      return rotation_axis_;
    }

    /** Get the fixed rotation matrix */
    mat3 <double> get_fixed_rotation() const {
      return fixed_rotation_;
    }

    /** Get the setting matrix */
    mat3 <double> get_setting_rotation() const {
      return setting_rotation_;
    }

    /** Get the number of scan points */
    std::size_t get_num_scan_points() const {
      return setting_rotation_at_scan_points_.size();
    }

    /** Get the setting rotation at scan points */
    scitbx::af::shared< mat3<double> > get_setting_rotation_at_scan_points() const {
      return setting_rotation_at_scan_points_;
    }

    /** Get the setting rotation at the scan point */
    mat3<double> get_setting_rotation_at_scan_point(std::size_t index) const {
      DXTBX_ASSERT(index < setting_rotation_at_scan_points_.size());
      return setting_rotation_at_scan_points_[index];
    }

    /** Set the rotation axis */
    void set_rotation_axis(vec3 <double> rotation_axis) {
      DXTBX_ASSERT(rotation_axis.length() > 0);
      rotation_axis_ = setting_rotation_.inverse() * rotation_axis.normalize();
    }

    /** Set the rotation axis */
    void set_rotation_axis_datum(vec3 <double> rotation_axis) {
      DXTBX_ASSERT(rotation_axis.length() > 0);
      rotation_axis_ = rotation_axis.normalize();
    }

    /** Set the fixed rotation matrix */
    void set_fixed_rotation(mat3 <double> fixed_rotation) {
      fixed_rotation_ = fixed_rotation;
    }

    /** Set the setting rotation matrix */
    void set_setting_rotation(mat3 <double> setting_rotation) {
      setting_rotation_ = setting_rotation;
    }

    /** Set the setting rotation matrix at scan-points */
    void set_setting_rotation_at_scan_points(const scitbx::af::const_ref< mat3<double> > &S) {
      setting_rotation_at_scan_points_ = scitbx::af::shared< mat3<double> >(S.begin(), S.end());
    }

    /** Reset the scan points */
    void reset_scan_points() {
      setting_rotation_at_scan_points_.clear();
    }

    /** Check rotation axes are (almost) the same */
    bool operator==(const Goniometer &rhs) const {
      double eps = 1.0e-6;

      // scan-varying model checks
      if (get_num_scan_points() > 0) {
        if (get_num_scan_points() != rhs.get_num_scan_points()) {
          return false;
        }
        for (std::size_t j = 0; j < get_num_scan_points(); ++j) {
          mat3<double> this_S = get_setting_rotation_at_scan_point(j);
          mat3<double> other_S = rhs.get_setting_rotation_at_scan_point(j);
          double d_S = 0.0;
          for (std::size_t i = 0; i < 9; ++i) {
            d_S += std::abs(this_S[i] - other_S[i]);
          }
          if (d_S > eps) {
            return false;
          }
        }
      }

      // static model checks
      return std::abs(angle_safe(rotation_axis_, rhs.rotation_axis_)) <= eps
      && fixed_rotation_.const_ref().all_approx_equal(rhs.fixed_rotation_.const_ref(), eps)
      && setting_rotation_.const_ref().all_approx_equal(rhs.setting_rotation_.const_ref(), eps);
    }

    /** Check rotation axes are not (almost) the same */
    bool operator!=(const Goniometer &goniometer) const {
      return !(*this == goniometer);
    }

    bool is_similar_to(const Goniometer &rhs,
                       double rotation_axis_tolerance,
                       double fixed_rotation_tolerance,
                       double setting_rotation_tolerance) const {

      // scan-varying model checks
      if (get_num_scan_points() != rhs.get_num_scan_points()) {
        return false;
      }
      for (std::size_t i = 0; i < get_num_scan_points(); ++i) {
        mat3<double> S_a = get_setting_rotation_at_scan_point(i);
        mat3<double> S_b = rhs.get_setting_rotation_at_scan_point(i);
        if (!S_a.const_ref().all_approx_equal(S_b.const_ref(), setting_rotation_tolerance)) {
          return false;
        }
      }

      // static model checks
      return std::abs(angle_safe(rotation_axis_, rhs.rotation_axis_)) <= rotation_axis_tolerance
      && fixed_rotation_.const_ref().all_approx_equal(rhs.fixed_rotation_.const_ref(), fixed_rotation_tolerance)
      && setting_rotation_.const_ref().all_approx_equal(rhs.setting_rotation_.const_ref(), setting_rotation_tolerance);
    }

    /** Rotate the goniometer about an axis */
    void rotate_around_origin(vec3<double> axis, double angle) {
      mat3<double> R = axis_and_angle_as_matrix(axis, angle);
      rotation_axis_ = rotation_axis_.rotate_around_origin(axis, angle);
      fixed_rotation_ = R * fixed_rotation_;
      setting_rotation_ = R * setting_rotation_;
    }

    friend std::ostream& operator<<(std::ostream& os, const Goniometer &gonio);

  protected:

    vec3 <double> rotation_axis_;
    mat3 <double> fixed_rotation_;
    mat3 <double> setting_rotation_;
    scitbx::af::shared< mat3<double> > setting_rotation_at_scan_points_;
  };

  /** Print goniometer data */
  inline
  std::ostream& operator<<(std::ostream& os, const Goniometer &g) {
    os << "Goniometer:\n";
    os << "    Rotation axis:   " << g.get_rotation_axis().const_ref() << "\n";
    os << "    Fixed rotation:  " << g.get_fixed_rotation().const_ref() << "\n";
    os << "    Setting rotation:" << g.get_setting_rotation().const_ref() << "\n";
    return os;
  }

}} // namespace dxtbx::model

#endif // DXTBX_MODEL_GONIOMETER_H
