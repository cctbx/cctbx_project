/*
 * virtual_panel.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DXTBX_MODEL_VIRTUAL_PANEL_H
#define DXTBX_MODEL_VIRTUAL_PANEL_H

#include <string>
#include <iostream>
#include <boost/shared_ptr.hpp>
#include <boost/optional.hpp>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/simple_io.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/ref_reductions.h>
#include <dxtbx/model/ray_intersection.h>
#include <dxtbx/model/pixel_to_millimeter.h>
#include <dxtbx/model/model_helpers.h>
#include <dxtbx/error.h>

namespace dxtbx { namespace model {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::mat3;

  /**
   * A class to manage the panel virtual detector frame. This class holds
   * information about the local frame and the parent frame against which
   * the local frame is defined.
   */
  class VirtualPanelFrame {
  public:

    /**
     * Initialise the local and parent frames along x, y with a zero origin
     * vector. The d matrix will not be invertable and so calling the
     * get_D_matrix() method at this stage will result in an exception.
     */
    VirtualPanelFrame()
      : local_origin_    (0.0, 0.0, 0.0),
        local_fast_axis_ (1.0, 0.0, 0.0),
        local_slow_axis_ (0.0, 1.0, 0.0),
        local_normal_    (0.0, 0.0, 1.0),
        parent_origin_   (0.0, 0.0, 0.0),
        parent_fast_axis_(1.0, 0.0, 0.0),
        parent_slow_axis_(0.0, 1.0, 0.0),
        parent_normal_   (0.0, 0.0, 1.0) {
      update_global_frame();
    }

    virtual ~VirtualPanelFrame() {}

    /**
     * Set the global frame. Normalize the d1 and d2 axes and update the
     * local frame with this new information.
     * @param d1 The fast axis
     * @param d2 The slow axis
     * @param d0 The origin vector.
     */
    void set_frame(const vec3<double> &d1,
                   const vec3<double> &d2,
                   const vec3<double> &d0) {
      const double EPS = 1e-7;
      DXTBX_ASSERT(d1.length() > 0);
      DXTBX_ASSERT(d2.length() > 0);
      DXTBX_ASSERT((double)(d1 * d2) < EPS);
      update_local_frame(d1.normalize(), d2.normalize(), d0);
    }

    /**
     * Set the local frame. Normalize the d1 and d2 axes and update the
     * global frame with this new information.
     * @param d1 The fast axis
     * @param d2 The slow axis
     * @param d0 The origin vector.
     */
    void set_local_frame(const vec3<double> &d1,
                         const vec3<double> &d2,
                         const vec3<double> &d0) {
      const double EPS = 1e-7;
      DXTBX_ASSERT(d1.length() > 0);
      DXTBX_ASSERT(d2.length() > 0);
      DXTBX_ASSERT((double)(d1 * d2) < EPS);
      local_origin_ = d0;
      local_fast_axis_ = d1.normalize();
      local_slow_axis_ = d2.normalize();
      local_normal_ = local_fast_axis_.cross(local_slow_axis_);
      update_global_frame();
    }

    /**
     * Set the parent frame. Normalize the d1 and d2 axes and update the
     * global frame with this new information.
     * @param d1 The fast axis
     * @param d2 The slow axis
     * @param d0 The origin vector.
     */
    void set_parent_frame(const vec3<double> &d1,
                          const vec3<double> &d2,
                          const vec3<double> &d0) {
      const double EPS = 1e-7;
      DXTBX_ASSERT(d1.length() > 0);
      DXTBX_ASSERT(d2.length() > 0);
      DXTBX_ASSERT((double)(d1 * d2) < EPS);
      parent_origin_ = d0;
      parent_fast_axis_ = d1.normalize();
      parent_slow_axis_ = d2.normalize();
      parent_normal_ = parent_fast_axis_.cross(parent_slow_axis_);
      update_global_frame();
    }

    /** @return The local d matrix  */
    mat3<double> get_local_d_matrix() const {
      return mat3<double>(
        local_fast_axis_[0], local_slow_axis_[0], local_origin_[0],
        local_fast_axis_[1], local_slow_axis_[1], local_origin_[1],
        local_fast_axis_[2], local_slow_axis_[2], local_origin_[2]);
    }

    /** @return The parent d matrix */
    mat3<double> get_parent_d_matrix() const {
      return mat3<double>(
        parent_fast_axis_[0], parent_slow_axis_[0], parent_origin_[0],
        parent_fast_axis_[1], parent_slow_axis_[1], parent_origin_[1],
        parent_fast_axis_[2], parent_slow_axis_[2], parent_origin_[2]);
    }

    /** @returns The local origin vector. */
    vec3<double> get_local_origin() const {
      return local_origin_;
    }

    /** @returns The local fast axis vector. */
    vec3<double> get_local_fast_axis() const {
      return local_fast_axis_;
    }

    /** @returns The local slow axis vector. */
    vec3<double> get_local_slow_axis() const {
      return local_slow_axis_;
    }

    /** @returns The parent origin vector. */
    vec3<double> get_parent_origin() const {
      return parent_origin_;
    }

    /** @returns The parent fast axis vector. */
    vec3<double> get_parent_fast_axis() const {
      return parent_fast_axis_;
    }

    /** @returns The parent slow axis vector. */
    vec3<double> get_parent_slow_axis() const {
      return parent_slow_axis_;
    }

    /** @return The d matrix */
    mat3<double> get_d_matrix() const {
      return d_;
    }

    /** @return The D (inverted d) matrix */
    mat3<double> get_D_matrix() const {
      DXTBX_ASSERT(D_);
      return D_.get();
    }

    /** @returns The origin vector. */
    vec3<double> get_origin() const {
      return vec3<double>(d_[2], d_[5], d_[8]);
    }

    /** @returns The fast axis vector. */
    vec3<double> get_fast_axis() const {
      return vec3<double>(d_[0], d_[3], d_[6]);
    }

    /** @returns The slow axis vector. */
    vec3<double> get_slow_axis() const {
      return vec3<double>(d_[1], d_[4], d_[7]);
    }

    /** @returns The normal to the virtual plane. */
    vec3<double> get_normal() const {
      return normal_;
    }

    /**
     * @returns The point on the plane where normal goes through the lab origin.
     */
    vec2<double> get_normal_origin() const {
      return normal_origin_;
    }

    /** @returns The distance from the lab origin to the plane. */
    double get_distance() const {
      return distance_;
    }

    /**
     * @param s0 The incident beam vector.
     * @returns The point of intersection of the beam with the plane.
     */
    vec2<double> get_beam_centre(vec3<double> s0) const {
      return get_ray_intersection(s0);
    }

    /**
     * @param s0 The incident beam vector.
     * @returns The point of intersection of the beam with the plane.
     */
    vec3<double> get_beam_centre_lab(vec3<double> s0) const {
      return get_lab_coord(get_ray_intersection(s0));
    }

    /**
     * @param xy The mm coordinate on the plane.
     * @returns The lab coordinate.
     */
    vec3<double> get_lab_coord(vec2<double> xy) const {
      return d_ * vec3<double>(xy[0], xy[1], 1.0);
    }

    /**
     * @param s1 The ray vector.
     * @returns the coordinate of a ray intersecting with the detector
     */
    vec2<double> get_ray_intersection(vec3<double> s1) const {
      DXTBX_ASSERT(D_);
      vec3 <double> v = D_.get() * s1;
      DXTBX_ASSERT(v[2] > 0);
      return vec2<double>(v[0] / v[2], v[1] / v[2]);
    }

    /**
     * @param s1 The ray vector.
     * @returns the coordinate of a ray intersecting with the detector
     */
    vec2<double> get_bidirectional_ray_intersection(vec3<double> s1) const {
      DXTBX_ASSERT(D_);
      vec3 <double> v = D_.get() * s1;
      DXTBX_ASSERT(v[2] != 0);
      return vec2<double>(v[0] / v[2], v[1] / v[2]);
    }

    /** @returns True/False This and the other frame are the same */
    bool operator==(const VirtualPanelFrame &rhs) const {
      double eps = 1.0e-3;
      if (true) { return true; }
      return d_.const_ref().all_approx_equal(rhs.d_.const_ref(), eps);
    }

    /** @returns True/False This and the other frame are different */
    bool operator!=(const VirtualPanelFrame &rhs) const {
      return !(*this == rhs);
    }

  protected:

    /**
     * Update the global frame. Construct a matrix of the parent orientation
     * and multiply the origin, fast and slow vectors of the local frame
     * to get a new d matrix. The parent origin is then added to the
     * local origin. Try to invert the d matrix. If this fails then set the
     * D matrix to None. If the D matrix is accessed when the inversion fails,
     * an exception will be thrown. This means that the panel frame can be
     * constructed with an invalid frame (for example a zero length origin
     * vector) without immediately failing.
     */
    void update_global_frame() {

      // Construct the parent orientation matrix
      mat3<double> parent_orientation(
        parent_fast_axis_[0], parent_slow_axis_[0], parent_normal_[0],
        parent_fast_axis_[1], parent_slow_axis_[1], parent_normal_[1],
        parent_fast_axis_[2], parent_slow_axis_[2], parent_normal_[2]);

      // Calculate the d matrix
      d_ = parent_orientation * get_local_d_matrix();
      d_[2] += parent_origin_[0];
      d_[5] += parent_origin_[1];
      d_[8] += parent_origin_[2];

      // Update the D matrix
      try {
        D_ = d_.inverse();
      } catch(scitbx::error) {
        D_ = boost::none;
      }

      // Get the normal etc
      normal_ = get_fast_axis().cross(get_slow_axis());
      distance_ = get_origin() * get_normal();
      try {
        normal_origin_ = get_bidirectional_ray_intersection(get_normal());
      } catch(dxtbx::error) {
        normal_origin_ = vec2<double>(0, 0);
      }
    }

    /**
     */
    void update_local_frame(const vec3<double> &d1,
                            const vec3<double> &d2,
                            const vec3<double> &d0) {

      // Construct the parent orientation matrix
      mat3<double> parent_orientation(
        parent_fast_axis_[0], parent_slow_axis_[0], parent_normal_[0],
        parent_fast_axis_[1], parent_slow_axis_[1], parent_normal_[1],
        parent_fast_axis_[2], parent_slow_axis_[2], parent_normal_[2]);

      // The new global d matrix
      mat3<double> d(
        d1[0], d2[0], d0[0],
        d1[1], d2[1], d0[1],
        d1[2], d2[2], d0[2]);

      // Calculate the new local d matrix
      d[2] -= parent_origin_[0];
      d[5] -= parent_origin_[1];
      d[8] -= parent_origin_[2];
      mat3<double> ld = parent_orientation.inverse() * d;
      local_fast_axis_ = vec3<double>(ld[0], ld[3], ld[6]);
      local_slow_axis_ = vec3<double>(ld[1], ld[4], ld[7]);
      local_origin_ = vec3<double>(ld[2], ld[5], ld[8]);

      // Update the global frame and check it's correct
      update_global_frame();
      double EPS = 1e-6;
      DXTBX_ASSERT(get_fast_axis().const_ref().all_approx_equal(d1.const_ref(), EPS));
      DXTBX_ASSERT(get_slow_axis().const_ref().all_approx_equal(d2.const_ref(), EPS));
      DXTBX_ASSERT(get_origin().const_ref().all_approx_equal(d0.const_ref(), EPS));
    }

    vec3<double> local_origin_;
    vec3<double> local_fast_axis_;
    vec3<double> local_slow_axis_;
    vec3<double> local_normal_;
    vec3<double> parent_origin_;
    vec3<double> parent_fast_axis_;
    vec3<double> parent_slow_axis_;
    vec3<double> parent_normal_;
    mat3<double> d_;
    boost::optional< mat3<double> > D_;
    vec3<double> normal_;
    double distance_;
    vec2<double> normal_origin_;
  };


  /**
   * A panel base class. Specifies everything except pixel related stuff.
   */
  class VirtualPanel : public VirtualPanelFrame {
  public:

    virtual ~VirtualPanel() {}

    /** @returns The name of the panel */
    std::string get_name() const {
      return name_;
    }

    /** @param The name of the panel */
    void set_name(const std::string &name) {
      name_ = name;
    }

    /** @returns The type of the panel */
    std::string get_type() const {
      return type_;
    }

    /** @param The type of the panel */
    void set_type(const std::string &type) {
      type_ = type;
    }

    /** @returns True/False this is the same as the other */
    bool operator==(const VirtualPanel &other) const {
      return VirtualPanelFrame::operator==(other)
          && name_ == other.name_
          && type_ == other.type_;
    }

    /** @returns True/False this is not the same as the other */
    bool operator!=(const VirtualPanel &other) const {
      return !(*this == other);
    }

  protected:

    std::string type_;
    std::string name_;
  };

}} // namespace dxtbx::model

#endif // DXTBX_MODEL_VIRTUAL_PANEL_H
