/*
* detector.h
*
*  Copyright (C) 2013 Diamond Light Source
*
*  Author: James Parkhurst
*
*  This code is distributed under the BSD license, a copy of which is
*  included in the root directory of this package.
*/
#ifndef DXTBX_MODEL_DETECTOR_H
#define DXTBX_MODEL_DETECTOR_H

#include <string>
//#include <boost/geometry.hpp>
//#include <boost/geometry/geometries/point.hpp>
//#include <boost/geometry/geometries/polygon.hpp>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/shared.h>
#include <dxtbx/error.h>
#include "panel.h"

namespace dxtbx { namespace model {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::mat3;
  using scitbx::af::double6;
  using scitbx::af::int4;

  // int4 array type
  typedef scitbx::af::flex<int4>::type flex_int4;
  typedef scitbx::af::shared<int4> shared_int4;
  typedef scitbx::af::flex<mat3<double> >::type flex_mat3_double;
  typedef scitbx::af::flex<std::string>::type flex_string;

  /**
  * A class representing a detector made up of multiple flat panel detectors.
  * The detector elements can be accessed in the same way as an array:
  *  detector[0] -> 1st detector panel.
  */
  class DetectorBase {
  public:

    typedef std::pair<int, vec2<double> > coord_type;

    // Panel list typedefs
    typedef Panel panel_type;
    typedef scitbx::af::shared <panel_type> panel_list_type;
    typedef panel_list_type::iterator iterator;

    /** Default constructor */
    DetectorBase() {}

    /** Initialise the detector with singla panel */
    DetectorBase(const Panel &panel) {
      panel_list_.push_back(panel);
    }

    /** Initialise the detector with an array of panels */
    DetectorBase(const panel_list_type &panel_list) {
      for (std::size_t i = 0; i < panel_list.size(); ++i) {
        panel_list_.push_back(panel_list[i]);
      }
    }

    /** Virtual destructor */
    virtual ~DetectorBase() {}

    /** Get the begin iterator */
    iterator begin() {
      return panel_list_.begin();
    }

    /** Get the end iterator */
    iterator end() {
      return panel_list_.end();
    }

    /** Add a panel to the list of panels */
    void add_panel(const panel_type &panel) {
      panel_list_.push_back(panel);
    }

    /** Remove all the panels */
    void remove_panels() {
      panel_list_.erase(panel_list_.begin(), panel_list_.end());
    }

    /** Remove a single panel */
    void remove_panel(std::size_t i) {
      panel_list_.erase(panel_list_.begin() + i);
    }

    /** Get the number of panels */
    std::size_t num_panels() const {
      return panel_list_.size();
    }

    /** Return a reference to a panel */
    panel_type& operator[](std::size_t index) {
      return panel_list_[index];
    }

    /** Return a const reference to a panel */
    const panel_type& operator[](std::size_t index) const {
      return panel_list_[index];
    }

    /** Check the detector panels are the same */
    bool operator==(const DetectorBase &detector) const {
      bool same = panel_list_.size() == detector.panel_list_.size();
      if (same) {
        for (std::size_t i = 0; i < panel_list_.size(); ++i) {
          same = same && (panel_list_[i] == detector.panel_list_[i]);
        }
      }
      return same;
    }

    /** Check the detector panels are not the same */
    bool operator!=(const DetectorBase &detector) const {
      return !(*this == detector);
    }

    /** Get a list of names */
    flex_string get_names() const {
      flex_string s(panel_list_.size());
      for (std::size_t i = 0; i < panel_list_.size(); ++i) {
        s[i] = panel_list_[i].get_name();
      }
      return s;
    }

    /** Get an array of d matrices from the panel list */
    flex_mat3_double get_d_matrices() const {
      flex_mat3_double d(panel_list_.size());
      for (std::size_t i = 0; i < panel_list_.size(); ++i) {
        d[i] = panel_list_[i].get_d_matrix();
      }
      return d;
    }

    /** Get an array of D matrices from the panel list */
    flex_mat3_double get_D_matrices() const {
      flex_mat3_double D(panel_list_.size());
      for (std::size_t i = 0; i < panel_list_.size(); ++i) {
        D[i] = panel_list_[i].get_D_matrix();
      }
      return D;
    }

    /** Get the maximum resolution of the detector */
    double get_max_resolution(vec3<double> s0, double wavelength) const {
      double d_min = 0;
      for (std::size_t i = 0; i < panel_list_.size(); ++i) {
        double d = panel_list_[i].get_max_resolution_at_corners(s0, wavelength);
        if (d < d_min || i == 0) d_min = d;
      }
      return d_min;
    }

    /** Get ray intersection with detector */
    coord_type get_ray_intersection(vec3<double> s1) const {
      coord_type pxy(-1, vec2<double>(0, 0));
      double w_max = 0;

      // Loop through all detectors. If the w component of the (u, v, w)
      // vector points in the correct direction, then calculate the coordinate.
      // If the coordinate is valid and the w component is greater than that of
      // the current closest valid coordinate, then set this coordinate as the
      // current best bet.
      for (std::size_t i = 0; i < panel_list_.size(); ++i) {
        vec3 <double> v = panel_list_[i].get_D_matrix() * s1;
        if (v[2] > 0) {
          vec2<double> xy_temp(v[0] / v[2], v[1] / v[2]);
          if (panel_list_[i].is_coord_valid_mm(xy_temp) && v[2] > w_max) {
            pxy = coord_type(i, xy_temp);
            w_max = v[2];
          }
        }
      }

      // If no coordinate was found then raise an exception
      // otherwise return the coordinate.
      DXTBX_ASSERT(w_max > 0);
      return pxy;
    }

//    /** Check if any panels intersect */
//    bool do_panels_intersect() const {
//      for (std::size_t j = 0; j < panel_list_.size()-1; ++j) {
//        for (std::size_t i = j+1; i < panel_list_.size(); ++i) {
//          if (panels_intersect(panel_list_[j], panel_list_[i])) {
//            return true;
//          }
//        }
//      }
//      return false;
//    }

    friend std::ostream& operator<< (std::ostream& , const DetectorBase&);

  protected:

//    /**
//     * Check if the detector planes intersect.
//     * @param a The first detector
//     * @param b The second detector
//     * @returns True/False do the detector planes intersect?
//     */
//    static bool
//    panels_intersect(const Panel &a, const Panel &b) {

//      using namespace boost::geometry;

//      typedef boost::geometry::model::point <double, 3, cs::cartesian> point;
//      typedef boost::geometry::model::polygon <point> polygon;

//      // Get the rectange of detector points
//      double6 rect_a(0, 0, 0, 0, 0, 0);// = a.get_image_rectangle();
//      double6 rect_b(0, 0, 0, 0, 0, 0);// = b.get_image_rectangle();

//      // Create a polygon for the panel a plane
//      polygon poly_a;
//      append(poly_a, point(rect_a[0], rect_a[1], rect_a[2]));
//      append(poly_a, point(rect_a[3], rect_a[1], rect_a[5]));
//      append(poly_a, point(rect_a[3], rect_a[4], rect_a[5]));
//      append(poly_a, point(rect_a[0], rect_a[4], rect_a[2]));
//      append(poly_a, point(rect_a[0], rect_a[1], rect_a[2]));

//      // Create a polygon for the panel b plane
//      polygon poly_b;
//      append(poly_b, point(rect_b[0], rect_b[1], rect_b[2]));
//      append(poly_b, point(rect_b[3], rect_b[1], rect_b[5]));
//      append(poly_b, point(rect_b[3], rect_b[4], rect_b[5]));
//      append(poly_b, point(rect_b[0], rect_b[4], rect_b[2]));
//      append(poly_b, point(rect_b[0], rect_b[1], rect_b[2]));

//      // Check if the polygons intersect
//      return intersects(poly_a, poly_b);
//    }

    panel_list_type panel_list_;
  };

  /** Print the detector information to the ostream */
  inline
  std::ostream& operator<< (std::ostream &os, const DetectorBase &d) {
    os << "Detector:\n";
    for (std::size_t i = 0; i < d.num_panels(); ++i) {
      std::stringstream ss;
      ss << d[i];
      std::string panel;
      std::string line;
      while (std::getline(ss, line)) {
        panel += "    " + line + "\n";
      }
      os << panel;
    }
    return os;
  }

  /** Add interface to detector that is valid if only 1 panel is present */
  class Detector : public DetectorBase {
  public:

    /** Default constructor */
    Detector() {}

    /** Initialise the detector with singla panel */
    Detector(const Panel &panel)
      : DetectorBase(panel) {}

    /** Initialise the detector with an array of panels */
    Detector(const panel_list_type &panel_list)
      : DetectorBase(panel_list) {}

    /** Virtual destructor */
    virtual ~Detector() {}

    /** Get the first panel type */
    std::string get_type() const {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      return DetectorBase::panel_list_[0].get_type();
    }

    /** Set the first panel type */
    void set_type(std::string type) {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      DetectorBase::panel_list_[0].set_type(type);
    }

    /** Get the first panel name */
    std::string get_name() const {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      return DetectorBase::panel_list_[0].get_name();
    }

    /** Set the first panel name */
    void set_name(std::string name) {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      DetectorBase::panel_list_[0].set_name(name);
    }

    /** Get the first panel fast axis */
    vec3 <double> get_fast_axis() const {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      return DetectorBase::panel_list_[0].get_fast_axis();
    }

    /** Get the first panel slow axis */
    vec3 <double> get_slow_axis() const {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      return DetectorBase::panel_list_[0].get_slow_axis();
    }

    /** Get the first panel origin */
    vec3 <double> get_origin() const {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      return DetectorBase::panel_list_[0].get_origin();
    }

    /** Get the first panel normal */
    vec3 <double> get_normal() const {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      return DetectorBase::panel_list_[0].get_normal();
    }

    /** Get the first panel pixel size */
    vec2 <double> get_pixel_size() const {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      return DetectorBase::panel_list_[0].get_pixel_size();
    }

    /** Set the first panel pixel size */
    void set_pixel_size(vec2 <double> pixel_size) {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      DetectorBase::panel_list_[0].set_pixel_size(pixel_size);
    }

    /** Get the first panel image size */
    vec2 <std::size_t> get_image_size() const {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      return DetectorBase::panel_list_[0].get_image_size();
    }

    /** Set the first panel image size */
    void set_image_size(vec2 <std::size_t> image_size) {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      DetectorBase::panel_list_[0].set_image_size(image_size);
    }

    /** Get the first panel trusted range */
    vec2 <double> get_trusted_range() const {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      return DetectorBase::panel_list_[0].get_trusted_range();
    }

    /** Set the first panel trusted range */
    void set_trusted_range(vec2 <double> trusted_range) {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      DetectorBase::panel_list_[0].set_trusted_range(trusted_range);
    }

    /** Get the mask array */
    shared_int4 get_mask() const {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      return DetectorBase::panel_list_[0].get_mask();
    }

    /** Set the first panel mask */
    void set_mask(const shared_int4 &mask) {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      DetectorBase::panel_list_[0].set_mask(mask);
    }

    /** Add an element to the first panel mask */
    void add_mask(int f0, int s0, int f1, int s1) {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      DetectorBase::panel_list_[0].add_mask(f0, s0, f1, s1);
    }

    /** Get the matrix of the first panel  coordinate system */
    mat3 <double> get_d_matrix() const {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      return DetectorBase::panel_list_[0].get_d_matrix();
    }

    /** Get the first panel inverse d matrix */
    mat3 <double> get_D_matrix() const {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      return DetectorBase::panel_list_[0].get_D_matrix();
    }

    /** Set the first panel origin, fast axis and slow axis */
    void set_frame(vec3<double> fast_axis, vec3<double> slow_axis,
        vec3<double> origin) {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      DetectorBase::panel_list_[0].set_frame(fast_axis, slow_axis, origin);
    }

    /** Get the distance from the sample to the first panel plane */
    double get_distance() const {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      return DetectorBase::panel_list_[0].get_distance();
    }

    /** Get the origin of the normal on the detector */
    vec2<double> get_normal_origin() const {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      return DetectorBase::panel_list_[0].get_normal_origin();
    }

    /** Get the first panel image size in millimeters */
    vec2<double> get_image_size_mm() const {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      return DetectorBase::panel_list_[0].get_image_size_mm();
    }

    /** Check the first panel value is valid */
    bool is_value_in_trusted_range(double value) const {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      return DetectorBase::panel_list_[0].is_value_in_trusted_range(value);
    }

    /** Check the first panel coordinate is valid */
    bool is_coord_valid(vec2<double> xy) const {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      return DetectorBase::panel_list_[0].is_coord_valid(xy);
    }

    /** Check the first panel coordinate is valid */
    bool is_coord_valid_mm(vec2<double> xy) const {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      return DetectorBase::panel_list_[0].is_coord_valid_mm(xy);
    }

    /** Get the beam centre in mm in the first panel  basis */
    vec2<double> get_beam_centre(vec3<double> s0) const {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      return DetectorBase::panel_list_[0].get_beam_centre(s0);
    }

    /** Get the first panel beam centre in lab coordinates */
    vec3<double> get_beam_centre_lab(vec3<double> s0) const {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      return DetectorBase::panel_list_[0].get_beam_centre_lab(s0);
    }

    /** Get the resolution at a given first panel pixel. */
    double get_resolution_at_pixel(vec3<double> s0, double wavelength,
        vec2<double> xy) const {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      return DetectorBase::panel_list_[0].get_resolution_at_pixel(s0,
        wavelength, xy);
    }

    /** Get the maximum resolution of the first panel. */
    double
    get_max_resolution_at_corners(vec3<double> s0, double wavelength) const {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      return DetectorBase::panel_list_[0].get_max_resolution_at_corners(s0,
        wavelength);
    }

    /** Get the maximum resolution of a full circle on the first panel. */
    double get_max_resolution_elipse(vec3<double> s0, double wavelength) const {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      return DetectorBase::panel_list_[0].get_max_resolution_elipse(s0,
        wavelength);
    }

    /** Get the first panel point (in mm) in lab coordinates */
    vec3<double> get_lab_coord(vec2<double> xy) const {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      return DetectorBase::panel_list_[0].get_lab_coord(xy);
    }

    /** Get the first panel point (in mm) in lab coordinates */
    vec3<double> get_pixel_lab_coord(vec2<double> xy) const {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      return DetectorBase::panel_list_[0].get_pixel_lab_coord(xy);
    }

    /** Map first panel coordinates in mm to pixels */
    vec2<double> millimeter_to_pixel(vec2<double> xy) const {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      return DetectorBase::panel_list_[0].millimeter_to_pixel(xy);
    }

    /** Map the first panel coordinates in pixels to millimeters */
    vec2<double> pixel_to_millimeter(vec2<double> xy) const {
      DXTBX_ASSERT(DetectorBase::panel_list_.size() == 1);
      return DetectorBase::panel_list_[0].pixel_to_millimeter(xy);
    }
  };

}} // namespace dxtbx::model

#endif // DXTBX_MODEL_DETECTOR_H
