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
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/simple_io.h>
#include <scitbx/array_family/ref_reductions.h>
#include <dxtbx/model/pixel_to_millimeter.h>
#include <dxtbx/model/model_helpers.h>
#include <dxtbx/model/virtual_panel.h>
#include <dxtbx/model/panel_data.h>
#include <dxtbx/model/pixel_to_millimeter.h>
#include <dxtbx/error.h>

namespace dxtbx { namespace model {

  using scitbx::af::double4;
  using boost::shared_ptr;
  using scitbx::af::tiny;

  /**
   * A panel class.
   */
  class Panel : public PanelData {
  public:

    /** Construct the panel with the simple px->mm strategy */
    Panel()
      : convert_coord_(new SimplePxMmStrategy()) {}

    /** Construct with data but no px/mm strategy */
    Panel(std::string type,
          std::string name,
          tiny<double,3> fast_axis,
          tiny<double,3> slow_axis,
          tiny<double,3> origin,
          tiny<double,2> pixel_size,
          tiny<std::size_t,2> image_size,
          tiny<double,2> trusted_range)
      : PanelData(type, name,
          fast_axis, slow_axis, origin,
          pixel_size, image_size,
          trusted_range),
        convert_coord_(new SimplePxMmStrategy()) {}

    /** Construct with data with px/mm strategy */
    Panel(std::string type,
          std::string name,
          tiny<double,3> fast_axis,
          tiny<double,3> slow_axis,
          tiny<double,3> origin,
          tiny<double,2> pixel_size,
          tiny<std::size_t,2> image_size,
          tiny<double,2> trusted_range,
          shared_ptr<PxMmStrategy> convert_coord)
      : PanelData(type, name,
          fast_axis, slow_axis, origin,
          pixel_size, image_size,
          trusted_range),
        convert_coord_(convert_coord) {}

    virtual ~Panel() {}

    /** Get the pixel to millimetre strategy */
    shared_ptr<PxMmStrategy> get_px_mm_strategy() const {
      return convert_coord_;
    }

    /** Set the pixel to millimetre strategy */
    void set_px_mm_strategy(shared_ptr<PxMmStrategy> strategy) {
      convert_coord_ = strategy;
    }

    /** Get the image size in millimeters */
    tiny<double,2> get_image_size_mm() const {
      return pixel_to_millimeter(tiny<double,2>(image_size_[0], image_size_[1]));
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
      tiny<double,2> size = get_image_size_mm();
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
    vec3<double> get_pixel_lab_coord(tiny<double,2> xy) const {
      tiny<double,2> xy_mm = pixel_to_millimeter(xy);
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
      double sintheta = sin(0.5 * angle_safe(s0, xyz));
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

    friend std::ostream& operator<<(std::ostream &os, const Panel &p);

  protected:
    shared_ptr<PxMmStrategy> convert_coord_;
  };

  /** Print panel information */
  inline
  std::ostream& operator<<(std::ostream &os, const Panel &p) {
    os << "Panel:" << std::endl;
    os << "  pixel_size:" << p.get_pixel_size().const_ref() << std::endl;
    os << "  image_size: " << p.get_image_size().const_ref() << std::endl;
    os << "  trusted_range: " << p.get_trusted_range().const_ref() << std::endl;
    os << "  fast_axis: " << p.get_fast_axis().const_ref() << std::endl;
    os << "  slow_axis: " << p.get_slow_axis().const_ref() << std::endl;
    os << "  origin: " << p.get_origin().const_ref() << std::endl;
    return os;
  }

}} // namespace dxtbx::model

#endif // DXTBX_MODEL_PANEL_H
