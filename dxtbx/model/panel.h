/*
* panel.h
*
*  Copyright (C) 2013 Diamond Light Source
*
*  Author: James Parkhurst
*
*  This code is distributed under the BSD license, a copy of which is
*  included in the root directory of this package.
*/
#ifndef DXTBX_MODEL_PANEL_H
#define DXTBX_MODEL_PANEL_H

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
#include <dxtbx/model/pixel_to_millimeter.h>
#include <dxtbx/model/model_helpers.h>
#include <dxtbx/error.h>

namespace dxtbx { namespace model {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::mat3;
  using scitbx::af::int4;
  using scitbx::af::double4;
  using boost::shared_ptr;

  /** Get the coordinate of a ray intersecting with the detector */
  inline
  vec2<double> plane_ray_intersection(mat3<double> D, vec3<double> s1) {
    vec3 <double> v = D * s1;
    DXTBX_ASSERT(v[2] > 0);
    return vec2<double>(v[0] / v[2], v[1] / v[2]);
  }

  /** Get the coordinate of a ray intersecting with the detector */
  inline
  vec2<double> bidirectional_plane_ray_intersection(
      mat3<double> D, vec3<double> s1) {
    vec3 <double> v = D * s1;
    DXTBX_ASSERT(v[2] != 0);
    return vec2<double>(v[0] / v[2], v[1] / v[2]);
  }

  /** Get world coordinate of plane xy */
  inline
  vec3<double> plane_world_coordinate(mat3<double> d, vec2<double> xy) {
    return d * vec3<double>(xy[0], xy[1], 1.0);
  }


  /**
   * A class to manage the panel virtual detector frame. This class holds
   * information about the local frame and the parent frame against which
   * the local frame is defined.
   */
  class PanelFrame {
  public:

    /**
     * Initialise the local and parent frames along x, y with a zero origin
     * vector. The d matrix will not be invertable and so calling the
     * get_D_matrix() method at this stage will result in an exception.
     */
    PanelFrame()
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

    virtual ~PanelFrame() {}

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
    bool operator==(const PanelFrame &other) const {
      double eps = 1.0e-3;
      return get_fast_axis().const_ref().all_approx_equal(
              other.get_fast_axis().const_ref(), eps)
          && get_slow_axis().const_ref().all_approx_equal(
              other.get_slow_axis().const_ref(), eps)
          && get_origin().const_ref().all_approx_equal(
              other.get_origin().const_ref(), eps);
    }

    /** @returns True/False This and the other frame are different */
    bool operator!=(const PanelFrame &other) const {
      return *this != other;
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
  class PanelBase : public PanelFrame {
  public:

    virtual ~PanelBase() {}

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
    bool operator==(const PanelBase &other) const {
      return PanelFrame::operator==(other)
          && name_ == other.name_
          && type_ == other.type_;
    }

    /** @returns True/False this is not the same as the other */
    bool operator!=(const PanelBase &other) const {
      return *this != other;
    }

  protected:

    std::string type_;
    std::string name_;
  };


  /**
   * A panel class.
   */
  class Panel : public PanelBase {
  public:

    /** Construct the panel with the simple px->mm strategy */
    Panel()
      : pixel_size_(0.0, 0.0),
        image_size_(0, 0),
        trusted_range_(0.0, 0.0),
        convert_coord_(new SimplePxMmStrategy()) {}

    /** Construct with data but no px/mm strategy */
    Panel(std::string type,
          std::string name,
          vec3 <double> fast_axis,
          vec3 <double> slow_axis,
          vec3 <double> origin,
          vec2 <double> pixel_size,
          vec2 <std::size_t> image_size,
          vec2 <double> trusted_range)
      : pixel_size_(pixel_size),
        image_size_(image_size),
        trusted_range_(trusted_range),
        convert_coord_(new SimplePxMmStrategy()) {
      set_type(type);
      set_name(name);
      set_local_frame(fast_axis, slow_axis, origin);
    }

    /** Construct with data with px/mm strategy */
    Panel(std::string type,
          std::string name,
          vec3 <double> fast_axis,
          vec3 <double> slow_axis,
          vec3 <double> origin,
          vec2 <double> pixel_size,
          vec2 <std::size_t> image_size,
          vec2 <double> trusted_range,
          shared_ptr<PxMmStrategy> convert_coord)
      : pixel_size_(pixel_size),
        image_size_(image_size),
        trusted_range_(trusted_range),
        convert_coord_(convert_coord) {
      set_type(type);
      set_name(name);
      set_local_frame(fast_axis, slow_axis, origin);
    }

    virtual ~Panel() {}

    /** Get the pixel size */
    vec2 <double> get_pixel_size() const {
      return pixel_size_;
    }

    /** Set the pixel size */
    void set_pixel_size(vec2 <double> pixel_size) {
      pixel_size_ = pixel_size;
    }

    /** Get the image size */
    vec2 <std::size_t> get_image_size() const {
      return image_size_;
    }

    /** Set the image size */
    void set_image_size(vec2 <std::size_t> image_size) {
      image_size_ = image_size;
    }

    /** Get the trusted range */
    vec2 <double> get_trusted_range() const {
      return trusted_range_;
    }

    /** Set the trusted range */
    void set_trusted_range(vec2 <double> trusted_range) {
      trusted_range_ = trusted_range;
    }

    /** Get the pixel to millimetre strategy */
    shared_ptr<PxMmStrategy> get_px_mm_strategy() const {
      return convert_coord_;
    }

    /** Set the pixel to millimetre strategy */
    void set_px_mm_strategy(shared_ptr<PxMmStrategy> strategy) {
      convert_coord_ = strategy;
    }

    /** Get the mask array */
    scitbx::af::shared<int4> get_mask() const {
      scitbx::af::shared<int4> result;
      std::copy(mask_.begin(), mask_.end(), std::back_inserter(result));
      return result;
    }

    /** Set the mask */
    void set_mask(const scitbx::af::const_ref<int4> &mask) {
      mask_.resize(0);
      std::copy(mask.begin(), mask.end(), mask_.begin());
    }

    /** Add an element to the mask */
    void add_mask(int f0, int s0, int f1, int s1) {
      mask_.push_back(int4(f0, f1, s0, s1));
    }

    /** Get the image size in millimeters */
    vec2<double> get_image_size_mm() const {
      return pixel_to_millimeter(vec2<double>(image_size_[0], image_size_[1]));
    }

    /** Check the value is valid */
    bool is_value_in_trusted_range(double value) const {
      return (trusted_range_[0] <= value && value < trusted_range_[1]);
    }

    /** Check the coordinate is valid */
    bool is_coord_valid(vec2<double> xy) const {
      return (0 <= xy[0] && xy[0] < image_size_[0])
          && (0 <= xy[1] && xy[1] < image_size_[1]);
    }

    /** Check the coordinate is valid */
    bool is_coord_valid_mm(vec2<double> xy) const {
      vec2<double> size = get_image_size_mm();
      return (0 <= xy[0] && xy[0] < size[0])
          && (0 <= xy[1] && xy[1] < size[1]);
    }

    vec2<double> get_normal_origin_px() const {
      return get_bidirectional_ray_intersection_px(get_normal());
    }

    /** Get the beam centre in mm in the detector basis */
    vec2<double> get_beam_centre_px(vec3<double> s0) const {
      return get_ray_intersection_px(s0);
    }

    /** Get the detector point (in mm) in lab coordinates */
    vec3<double> get_pixel_lab_coord(vec2<double> xy) const {
      vec2<double> xy_mm = pixel_to_millimeter(xy);
      return get_lab_coord(xy_mm);
    }

    /** Get the ray intersection in pixel coordinates */
    vec2<double> get_ray_intersection_px(vec3<double> s1) const {
      return millimeter_to_pixel(get_ray_intersection(s1));
    }

    /** Get the ray intersection in pixel coordinates */
    vec2<double> get_bidirectional_ray_intersection_px(vec3<double> s1) const {
      return millimeter_to_pixel(get_bidirectional_ray_intersection(s1));
    }

    /** Map coordinates in mm to pixels */
    vec2<double> millimeter_to_pixel(vec2<double> xy) const {
      return convert_coord_->to_pixel(*this, xy);
    }

    /** Map the coordinates in pixels to millimeters */
    vec2<double> pixel_to_millimeter(vec2<double> xy) const {
      return convert_coord_->to_millimeter(*this, xy);
    }

    /**
     * Get the resolution at a given pixel.
     * @param s0 The incident beam vector
     * @param xy The pixel coordinate
     * @returns The resolution at that point.
     */
    double get_resolution_at_pixel(vec3<double> s0, vec2<double> xy) const {
      DXTBX_ASSERT(s0.length() > 0);
      vec3<double> xyz = get_pixel_lab_coord(xy);
      vec3<double> beam_centre = get_beam_centre_lab(s0);
      double sintheta = sin(0.5 * angle_safe(beam_centre, xyz));
      DXTBX_ASSERT(sintheta != 0);
      return 1.0 / (2.0 * s0.length() * sintheta);
    }

    /**
     * Get the maximum resolution of the detector (i.e. look at each corner
     * and find the maximum resolution.)
     * @param beam The beam parameters
     * @returns The maximum resolution at the detector corners.
     */
    double
    get_max_resolution_at_corners(vec3<double> s0) const {
      int fast = image_size_[0], slow = image_size_[1];
      return scitbx::af::min(double4(
        get_resolution_at_pixel(s0, vec2<double>(0, 0)),
        get_resolution_at_pixel(s0, vec2<double>(0, slow)),
        get_resolution_at_pixel(s0, vec2<double>(fast, 0)),
        get_resolution_at_pixel(s0, vec2<double>(fast, slow))));
    }

    /**
     * Get the maximum resolution of a full circle on the detector. Get the
     * beam centre in pixels. Then find the coordinates on the edges making
     * a cross-hair with the beam centre. Calculate the resolution at these
     * corners and choose the minimum angle.
     * @param beam The beam parameters
     * @returns The maximum resolution at the detector corners.
     */
    double
    get_max_resolution_ellipse(vec3<double> s0) const {
      int fast = image_size_[0], slow = image_size_[1];
      vec2<double> c = get_beam_centre_px(s0);
      return scitbx::af::max(double4(
        get_resolution_at_pixel(s0, vec2<double>(0, c[1])),
        get_resolution_at_pixel(s0, vec2<double>(fast, c[1])),
        get_resolution_at_pixel(s0, vec2<double>(c[0], 0)),
        get_resolution_at_pixel(s0, vec2<double>(c[0], slow))));
    }

    /** @returns True/False this is the same as the other */
    bool operator==(const Panel &other) const {
      return PanelBase::operator==(other)
          && image_size_ == other.image_size_
          && pixel_size_.const_ref().all_approx_equal(
              other.pixel_size_.const_ref(), 1e-6)
          && trusted_range_.const_ref().all_approx_equal(
              other.trusted_range_.const_ref(), 1e-6);
    }

    /** @returns True/False this is not the same as the other */
    bool operator!=(const Panel &other) const {
      return !(*this == other);
    }

    friend std::ostream& operator<<(std::ostream &os, const Panel &p);

  protected:
    vec2 <double> pixel_size_;
    vec2 <std::size_t> image_size_;
    vec2 <double> trusted_range_;

    scitbx::af::shared<int4> mask_;

    // Pixel to millimeter function
    shared_ptr<PxMmStrategy> convert_coord_;
  };

  /** Print panel information */
  inline
  std::ostream& operator<<(std::ostream &os, const Panel &p) {
    os << "Panel:\n";
    os << "  pixel_size: \n" << p.get_pixel_size().const_ref();
    os << "  image_size: \n" << p.get_image_size().const_ref();
    os << "  trusted_range: \n" << p.get_trusted_range().const_ref();
    os << "  fast_axis: \n" << p.get_fast_axis().const_ref();
    os << "  slow_axis: \n" << p.get_slow_axis().const_ref();
    os << "  origin: \n" << p.get_origin().const_ref();
    return os;
  }

}} // namespace dxtbx::model

#endif // DXTBX_MODEL_PANEL_H
