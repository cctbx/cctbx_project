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
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/array_family/simple_io.h>
#include <scitbx/array_family/simple_tiny_io.h>
#include <dxtbx/error.h>
#include <dxtbx/model/pixel_to_millimeter.h>
#include "model_helpers.h"

namespace dxtbx { namespace model {

  using boost::shared_ptr;
  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::mat3;
  using scitbx::af::int4;
  using scitbx::af::double4;
  using scitbx::af::double6;

  // int4 array type
  typedef scitbx::af::flex<int4>::type flex_int4;
  typedef scitbx::af::shared<int4> shared_int4;

  /**
   * A class representing a detector panel. A detector can have multiple
   * panels which are each represented by this class.
   */
  class Panel {
  public:

    /** The default constructor */
    Panel()
      : type_("Unknown"),
        d_(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        D_(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        pixel_size_(0.0, 0.0),
        image_size_(0, 0),
        trusted_range_(0, 0),
        normal_(0.0, 0.0, 0.0),
        normal_origin_(0.0, 0.0),
        distance_(0.0),
        convert_coord_(new SimplePxMmStrategy()) {}

    /**
    * Initialise the detector panel.
    * @param type The type of the detector panel
    * @param fast_axis The fast axis of the detector. The vector is normalized.
    * @param slow_axis The slow axis of the detector. The vector is normalized.
    * @param normal The detector normal. The given vector is normalized.
    * @param origin The detector origin
    * @param pixel_size The size of the individual pixels
    * @param image_size The size of the detector panel (in pixels)
    * @param trusted_range The trusted range of the detector pixel values.
    * @param distance The distance from the detector to the crystal origin
    */
    Panel(std::string type,
        vec3 <double> fast_axis,
        vec3 <double> slow_axis,
        vec3 <double> origin,
        vec2 <double> pixel_size,
        vec2 <std::size_t> image_size,
        vec2 <double> trusted_range)
      : convert_coord_(new SimplePxMmStrategy()) {
      DXTBX_ASSERT(fast_axis.length() > 0);
      DXTBX_ASSERT(slow_axis.length() > 0);
      type_ = type;
      d_ = create_d_matrix(fast_axis.normalize(),
        slow_axis.normalize(), origin);
      D_ = d_.inverse();
      pixel_size_ = pixel_size;
      image_size_ = image_size;
      trusted_range_ = trusted_range;

      // Update the normal etc
      update_normal();
    }

    /**
    * Initialise the detector panel.
    * @param type The type of the detector panel
    * @param fast_axis The fast axis of the detector. The vector is normalized.
    * @param slow_axis The slow axis of the detector. The vector is normalized.
    * @param normal The detector normal. The given vector is normalized.
    * @param origin The detector origin
    * @param pixel_size The size of the individual pixels
    * @param image_size The size of the detector panel (in pixels)
    * @param trusted_range The trusted range of the detector pixel values.
    * @param distance The distance from the detector to the crystal origin
    * @param px_mm The pixel to millimeter strategy
    */
    Panel(std::string type,
        vec3 <double> fast_axis,
        vec3 <double> slow_axis,
        vec3 <double> origin,
        vec2 <double> pixel_size,
        vec2 <std::size_t> image_size,
        vec2 <double> trusted_range,
        shared_ptr<PxMmStrategy> px_mm)
      : convert_coord_(px_mm) {
      DXTBX_ASSERT(fast_axis.length() > 0);
      DXTBX_ASSERT(slow_axis.length() > 0);
      type_ = type;
      d_ = create_d_matrix(fast_axis.normalize(),
        slow_axis.normalize(), origin);
      D_ = d_.inverse();
      pixel_size_ = pixel_size;
      image_size_ = image_size;
      trusted_range_ = trusted_range;

      // Update the normal etc
      update_normal();
    }

    /** Virtual destructor */
    virtual ~Panel() {}

    /** Get the sensor type */
    std::string get_type() const {
      return type_;
    }

    /** Set the detector panel type */
    void set_type(std::string type) {
      type_ = type;
    }

    /** Get the fast axis */
    vec3 <double> get_fast_axis() const {
      return vec3<double>(d_[0], d_[3], d_[6]);
    }

    /** Get the slow axis */
    vec3 <double> get_slow_axis() const {
      return vec3<double>(d_[1], d_[4], d_[7]);
    }

    /** Get the pixel origin */
    vec3 <double> get_origin() const {
      return vec3<double>(d_[2], d_[5], d_[8]);
    }

    /** Get the normal */
    vec3 <double> get_normal() const {
      return normal_;
    }

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

    /** Get the mask array */
    shared_int4 get_mask() const {
      return mask_;
    }

    /** Set the mask */
    void set_mask(const shared_int4 &mask) {
      mask_ = mask;
    }

    /** Add an element to the mask */
    void add_mask(int f0, int s0, int f1, int s1) {
      mask_.push_back(int4(f0, f1, s0, s1));
    }

    /** Get the matrix of the detector coordinate system */
    mat3 <double> get_d_matrix() const {
      return d_;
    }

    /** Get the inverse d matrix */
    mat3 <double> get_D_matrix() const {
      return D_;
    }

    /** Set the origin, fast axis and slow axis */
    void set_frame(vec3<double> fast_axis, vec3<double> slow_axis,
        vec3<double> origin) {
      DXTBX_ASSERT(fast_axis.length() > 0);
      DXTBX_ASSERT(slow_axis.length() > 0);
      d_ = create_d_matrix(
        fast_axis.normalize(),
        slow_axis.normalize(),
        origin);
      D_ = d_.inverse();

      // Update the normal etc
      update_normal();
    }

    /** Get the distance from the sample to the detector plane */
    double get_distance() const {
      return distance_;
    }

    /** Get the origin of the normal. */
    vec2<double> get_normal_origin() const {
      return normal_origin_;
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

    /** Get the beam centre in mm in the detector basis */
    vec2<double> get_beam_centre(vec3<double> s0) const {
      return get_ray_intersection(s0);
    }

    /** Get the beam centre in lab coordinates */
    vec3<double> get_beam_centre_lab(vec3<double> s0) const {
      return get_lab_coord(get_ray_intersection(s0));
    }

    /**
     * Get the resolution at a given pixel.
     * @param beam The beam parameters
     * @param xy The pixel coordinate
     * @returns The resolution at that point.
     */
    double get_resolution_at_pixel(vec3<double> s0, double wavelength,
        vec2<double> xy) const {

      // Get lab coordinates of detector corners
      vec3<double> xyz = get_pixel_lab_coord(xy);

      // Calculate the point at which the beam intersects with the detector
      vec3<double> beam_centre = get_beam_centre_lab(s0);

      // Calculate the angle to the point
      double sintheta = sin(0.5 * angle_safe(beam_centre, xyz));
      DXTBX_ASSERT(sintheta != 0);

      // Return d = lambda / (2sin(theta))
      return wavelength / (2.0 * sintheta);
    }

    /**
     * Get the maximum resolution of the detector (i.e. look at each corner
     * and find the maximum resolution.)
     * @param beam The beam parameters
     * @returns The maximum resolution at the detector corners.
     */
    double
    get_max_resolution_at_corners(vec3<double> s0, double wavelength) const {

      // Get lab coordinates of detector corners
      int fast = image_size_[0], slow = image_size_[1];
      vec3<double> xyz00 = get_origin();
      vec3<double> xyz01 = get_pixel_lab_coord(vec2<double>(0, slow));
      vec3<double> xyz10 = get_pixel_lab_coord(vec2<double>(fast, 0));
      vec3<double> xyz11 = get_pixel_lab_coord(vec2<double>(fast, slow));

      // Calculate the point at which the beam intersects with the detector
      vec3<double> beam_centre = get_beam_centre_lab(s0);

      // Calculate half the maximum angle to the corners
      double sintheta = sin(0.5 * scitbx::af::max(double4(
        angle_safe(beam_centre, xyz00), angle_safe(beam_centre, xyz01),
        angle_safe(beam_centre, xyz10), angle_safe(beam_centre, xyz11))));
      DXTBX_ASSERT(sintheta != 0);

      // Return d = lambda / (2sin(theta))
      return wavelength / (2.0 * sintheta);
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
    get_max_resolution_elipse(vec3<double> s0, double wavelength) const {

      // Get beam centre in pixels and get the coordinates with the centre
      // as a crosshair
      int fast = image_size_[0], slow = image_size_[1];
      vec2<double> c = millimeter_to_pixel(get_beam_centre(s0));
      vec3<double> xyz0c = get_pixel_lab_coord(vec2<double>(0.0, c[1]));
      vec3<double> xyz1c = get_pixel_lab_coord(vec2<double>(fast, c[1]));
      vec3<double> xyzc0 = get_pixel_lab_coord(vec2<double>(c[0], 0.0));
      vec3<double> xyzc1 = get_pixel_lab_coord(vec2<double>(c[0], slow));

      // Calculate the point at which the beam intersects with the detector
      vec3<double> beam_centre = get_beam_centre_lab(s0);

      // Calculate half the minimum angle to the sides around the beam centre
      double sintheta = sin(0.5 * scitbx::af::min(double4(
        angle_safe(beam_centre, xyz0c), angle_safe(beam_centre, xyzc0),
        angle_safe(beam_centre, xyz1c), angle_safe(beam_centre, xyzc1))));

      // Return d = lambda / (2sin(theta))
      return wavelength / (2.0 * sintheta);
    }

    /** Get the detector point (in mm) in lab coordinates */
    vec3<double> get_lab_coord(vec2<double> xy) const {
      return d_ * vec3<double>(xy[0], xy[1], 1.0);
    }

    /** Get the detector point (in mm) in lab coordinates */
    vec3<double> get_pixel_lab_coord(vec2<double> xy) const {
      vec2<double> xy_mm = pixel_to_millimeter(xy);
      return get_lab_coord(xy_mm);
    }

    /** Get the coordinate of a ray intersecting with the detector */
    vec2<double> get_ray_intersection(vec3<double> s1) const {
      vec3 <double> v = D_ * s1;
      DXTBX_ASSERT(v[2] > 0);
      return vec2<double>(v[0] / v[2], v[1] / v[2]);
    }

    /** Get the coordinate of a ray intersecting with the detector */
    vec2<double> get_bidirectional_ray_intersection(vec3<double> s1) const {
      vec3 <double> v = D_ * s1;
      DXTBX_ASSERT(v[2] != 0);
      return vec2<double>(v[0] / v[2], v[1] / v[2]);
    }

    /** Get the ray intersection in pixel coordinates */
    vec2<double> get_ray_intersection_px(vec3<double> s1) const {
      return millimeter_to_pixel(get_ray_intersection(s1));
    }

    /** Map coordinates in mm to pixels */
    vec2<double> millimeter_to_pixel(vec2<double> xy) const {
      return convert_coord_->to_pixel(*this, xy);
    }

    /** Map the coordinates in pixels to millimeters */
    vec2<double> pixel_to_millimeter(vec2<double> xy) const {
      return convert_coord_->to_millimeter(*this, xy);
    }

    /** Check the detector axis basis vectors are (almost) the same */
    bool operator==(const Panel &panel) const {
      double eps = 1.0e-3;
      return std::abs(angle_safe(get_fast_axis(), panel.get_fast_axis())) < eps
          && std::abs(angle_safe(get_slow_axis(), panel.get_slow_axis())) < eps
          && std::abs(angle_safe(get_origin(), panel.get_origin())) < eps
          && image_size_ == panel.image_size_;
    }

    /** Check the detector axis basis vectors are not (almost) the same */
    bool operator!=(const Panel &detector) const {
      return !(*this == detector);
    }

    friend std::ostream& operator<< (std::ostream& , const Panel&);

  protected:

    static mat3<double> create_d_matrix(vec3<double> fast_axis,
        vec3<double> slow_axis, vec3<double> origin) {
      return mat3 <double> (
        fast_axis[0], slow_axis[0], origin[0],
        fast_axis[1], slow_axis[1], origin[1],
        fast_axis[2], slow_axis[2], origin[2]);
    }

    /** Update the normal properties. */
    void update_normal() {
      normal_ = get_fast_axis().cross(get_slow_axis());
      normal_origin_ = get_bidirectional_ray_intersection(get_normal());
      distance_ = get_origin() * get_normal();
    }

    std::string type_;
    mat3<double> d_;
    mat3<double> D_;
    vec2 <double> pixel_size_;
    vec2 <std::size_t> image_size_;
    vec2 <double> trusted_range_;
    shared_int4 mask_;
    vec3<double> normal_;
    vec2<double> normal_origin_;
    double distance_;

    // Pixel to millimeter function
    shared_ptr<PxMmStrategy> convert_coord_;
  };

  /** Print the panel information to the ostream */
  inline
  std::ostream& operator<< (std::ostream &os, const Panel &p) {
    os << "Panel:\n";
    os << "    type:          " << p.get_type() << "\n";
    os << "    fast axis:     " << p.get_fast_axis().const_ref() << "\n";
    os << "    slow axis:     " << p.get_slow_axis().const_ref() << "\n";
    os << "    origin:        " << p.get_origin().const_ref() << "\n";
    os << "    normal:        " << p.get_normal().const_ref() << "\n";
    os << "    pixel size:    " << p.get_pixel_size().const_ref() << "\n";
    os << "    image size:    " << p.get_image_size().const_ref() << "\n";
    os << "    trusted range: " << p.get_trusted_range().const_ref() << "\n";
    return os;
  }

}} // namespace dxtbx::model

#endif // DXTBX_MODEL_PANEL_H
