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

namespace dxtbx { namespace model {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::mat3;
  using scitbx::af::double6;
  using scitbx::af::int4;

  // int4 array type
  typedef scitbx::af::flex<int4>::type flex_int4;
  typedef scitbx::af::shared<int4> shared_int4;

  /** A base class for detectors */
  class DetectorBase {};

  /**
  * A class representing a detector panel. A detector can have multiple
  * panels which are each represented by this class.
  *
  * The class contains the following members accessible through getter and
  * setter methods:
  *   type A string representing the type of the detector panel (i.e the
  *     manufactors name or some such identifier).
  *
  *   x_axis A unit vector pointing along the x (fast) axis of the panel. The
  *     vector is given with respect to the laboratory coordinate frame.
  *
  *   y_axis A unit vector pointing along the y (slow) axis of the panel. The
  *     vector is given with respect to the laboratory coordinate frame.
  *
  *   normal A unit vector giving the normal to the panel plane. The
  *     vector is given with respect to the laboratory coordinate frame.
  *
  *   origin The origin of the detector plane. (i.e. the laboratory coordinate
  *     of the edge of the zeroth pixel
  *
  *   pixel_size The size of the pixels in mm. The convention used is that
  *     of (y, x) i.e. (slow, fast)
  *
  *   image_size The size of the panel in pixels. The convention used is that
  *     of (y, x) i.e. (slow, fast)
  *
  *   trusted_range The range of counts that are considered reliable, given
  *     in the range [min, max].
  *
  *   distance The signed distance from the source to the detector
  *
  * In the document detailing the conventions used:
  *    type -> *unspecified*
  *    x_axis -> d1
  *    y_axis -> d2
  *    normal -> d3
  *    origin -> *unspecified*
  *    pixel_size -> *unspecified*
  *    image_size -> *unspecified*
  *    trusted_range -> *unspecified*
  *    distance -> *unspecified*
  */
  class Detector : public DetectorBase {
  public:

    /** The default constructor */
    Detector()
      : type_("Unknown"),
      fast_axis_(1.0, 0.0, 0.0),
      slow_axis_(0.0, 1.0, 0.0),
      normal_(0.0, 0.0, 1.0),
      origin_(0.0, 0.0, 0.0),
      pixel_size_(0.0, 0.0),
      image_size_(0, 0),
      trusted_range_(0, 0) {}

    /**
    * Initialise the detector panel.
    * @param type The type of the detector panel
    * @param fast_axis The fast axis of the detector. The given vector is normalized.
    * @param slow_axis The slow axis of the detector. The given vector is normalized.
    * @param normal The detector normal. The given vector is normalized.
    * @param origin The detector origin
    * @param pixel_size The size of the individual pixels
    * @param image_size The size of the detector panel (in pixels)
    * @param trusted_range The trusted range of the detector pixel values.
    * @param distance The distance from the detector to the crystal origin
    */
    Detector(std::string type,
        vec3 <double> fast_axis,
        vec3 <double> slow_axis,
        vec3 <double> origin,
        vec2 <double> pixel_size,
        vec2 <std::size_t> image_size,
        vec2 <double> trusted_range)
      : type_(type),
      fast_axis_(fast_axis.normalize()),
      slow_axis_(slow_axis.normalize()),
      normal_((fast_axis_.cross(slow_axis_)).normalize()),
      origin_(origin),
      pixel_size_(pixel_size),
      image_size_(image_size),
      trusted_range_(trusted_range) {}

    /** Virtual destructor */
    virtual ~Detector() {}

    /** Get the sensor type */
    std::string get_type() const {
      return type_;
    }

    /** Get the fast axis */
    vec3 <double> get_fast_axis() const {
      return fast_axis_;
    }

    /** Get the slow axis */
    vec3 <double> get_slow_axis() const {
      return slow_axis_;
    }

    /** Get the normal */
    vec3 <double> get_normal() const {
      return normal_;
    }

    /** Get the pixel origin */
    vec3 <double> get_origin() const {
      return origin_;
    }

    /** Get the pixel size */
    vec2 <double> get_pixel_size() const {
      return pixel_size_;
    }

    /** Get the image size */
    vec2 <std::size_t> get_image_size() const {
      return image_size_;
    }

    /** Get the trusted range */
    vec2 <double> get_trusted_range() const {
      return trusted_range_;
    }

    /** Get the matrix of the detector coordinate system */
    mat3 <double> get_d_matrix() const {
      return mat3 <double> (
        fast_axis_[0], slow_axis_[0], origin_[0],
        fast_axis_[1], slow_axis_[1], origin_[1],
        fast_axis_[2], slow_axis_[2], origin_[2]);
    }

    /** Get the inverse d matrix */
    mat3 <double> get_D_matrix() const {
      return get_d_matrix().inverse();
    }

    /** Get the mask array */
    shared_int4 get_mask() const {
      return mask_;
    }

    /** Set the detector panel type */
    void set_type(std::string type) {
      type_ = type;
    }

    /** Set the fast axis */
    void set_fast_axis(vec3 <double> fast_axis) {
      fast_axis_ = fast_axis;
      normal_ = fast_axis_.cross(slow_axis_);
    }

    /** Set the slow axis */
    void set_slow_axis(vec3 <double> slow_axis) {
      slow_axis_ = slow_axis;
      normal_ = fast_axis_.cross(slow_axis_);
    }

    /** Set the origin */
    void set_origin(vec3 <double> origin) {
      origin_ = origin;
    }

    /** Set the pixel size */
    void set_pixel_size(vec2 <double> pixel_size) {
      pixel_size_ = pixel_size;
    }

    /** Set the image size */
    void set_image_size(vec2 <std::size_t> image_size) {
      image_size_ = image_size;
    }

    /** Set the trusted range */
    void set_trusted_range(vec2 <double> trusted_range) {
      trusted_range_ = trusted_range;
    }

    /** Set the matrix of the detector coordinate system */
    void set_d_matrix(mat3 <double> d) {
      fast_axis_ = d.get_column(0);
      slow_axis_ = d.get_column(1);
      origin_ = d.get_column(2);
    }

    /** Set the inverse d matrix */
    void set_D_matrix(mat3 <double> D) {
      set_d_matrix(D.inverse());
    }

    /** Set the mask */
    void set_mask(const shared_int4 &mask) {
      mask_ = mask;
    }

    void add_mask(int f0, int s0, int f1, int s1) {
      mask_.push_back(int4(f0, f1, s0, s1));
    }

    /** Check the detector axis basis vectors are (almost) the same */
    bool operator==(const Detector &detector) const {
      double eps = 1.0e-6;
      double d_fast = fast_axis_.angle(detector.fast_axis_);
      double d_slow = slow_axis_.angle(detector.slow_axis_);
      double d_origin = origin_.angle(detector.origin_);
      double d_size_slow =
        std::abs((int)image_size_[0] - (int)detector.image_size_[0]);
      double d_size_fast =
        std::abs((int)image_size_[1] - (int)detector.image_size_[1]);
      return d_fast <= eps && d_slow <= eps && d_origin <= eps &&
             d_size_slow <= eps && d_size_fast <= eps;
    }

    /** Check the detector axis basis vectors are not (almost) the same */
    bool operator!=(const Detector &detector) const {
      return !(*this == detector);
    }

    /** Get the pixel in lab coordinates */
    template <typename CoordType>
    vec3<double> get_pixel_lab_coord(CoordType xy) const {
      vec2<double> xy_mm = pixel_to_millimeter(xy);
      return origin_ + fast_axis_ * xy_mm[0] + slow_axis_ * xy_mm[1];
    }

    /** Get the image size in millimeters */
    vec2<double> get_image_size_mm() const {
      vec2<double> xy = pixel_to_millimeter(image_size_);
      return vec2<double>(std::abs(xy[0]), std::abs(xy[1]));
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

    /** Map coordinates in mm to pixels */
    template <typename CoordType>
    vec2<double> millimeter_to_pixel(CoordType xy) const {
      return vec2<double> (xy[0] / pixel_size_[0], xy[1] / pixel_size_[1]);
    }

    /** Map the coordinates in pixels to millimeters */
    template <typename CoordType>
    vec2<double> pixel_to_millimeter(CoordType xy) const {
      return vec2<double> (xy[0] * pixel_size_[0], xy[1] * pixel_size_[1]);
    }

  protected:

    std::string type_;
    vec3 <double> fast_axis_;
    vec3 <double> slow_axis_;
    vec3 <double> normal_;
    vec3 <double> origin_;
    vec2 <double> pixel_size_;
    vec2 <std::size_t> image_size_;
    vec2 <double> trusted_range_;
    shared_int4 mask_;
  };

  /**
  * A class representing a detector made up of multiple flat panel detectors.
  * The detector elements can be accessed in the same way as an array:
  *  detector[0] -> 1st detector panel.
  */
  class MultiPanelDetector : public DetectorBase {
  public:

    typedef std::pair<int, vec2<double> > coordinate_type;

    // Panel list typedefs
    typedef Detector panel_type;
    typedef scitbx::af::shared <panel_type> panel_list_type;
    typedef panel_list_type::iterator iterator;

    /** Default constructor */
    MultiPanelDetector()
      : type_("Unknown") {}

    /** Initialise the detector */
    MultiPanelDetector(std::string type)
      :type_(type) {}

    /** Virtual destructor */
    virtual ~MultiPanelDetector() {}

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
    bool operator==(const MultiPanelDetector &detector) {
      bool same = panel_list_.size() == detector.panel_list_.size();
      if (same) {
        for (std::size_t i = 0; i < panel_list_.size(); ++i) {
          same = same && (panel_list_[i] == detector.panel_list_[i]);
        }
      }
      return same;
    }

    /** Check the detector panels are not the same */
    bool operator!=(const MultiPanelDetector &detector) {
      return !(*this == detector);
    }

    /** Check the value is valid */
    bool is_value_in_trusted_range(int panel, double value) const {
      return 0 <= panel && panel < panel_list_.size()
          && panel_list_[panel].is_value_in_trusted_range(value);
    }

    /** Check the coordinate is valid */
    bool is_coord_valid(coordinate_type pxy) const {
      return 0 <= pxy.first && pxy.first < panel_list_.size()
          && panel_list_[pxy.first].is_coord_valid(pxy.second);
    }

    /** Map coordinates in mm to pixels */
    vec2<double> millimeter_to_pixel(coordinate_type pxy) const {
      return panel_list_[pxy.first].millimeter_to_pixel(pxy.second);
    }

    /** Map the coordinates in pixels to millimeters */
    vec2<double> pixel_to_millimeter(coordinate_type pxy) const {
      return panel_list_[pxy.first].pixel_to_millimeter(pxy.second);
    }

    /** Check if any panels intersect */
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

  protected:

    /**
     * Check if the detector planes intersect.
     * @param a The first detector
     * @param b The second detector
     * @returns True/False do the detector planes intersect?
     */
//    static bool
//    panels_intersect(const Detector &a, const Detector &b) {

//      using namespace boost::geometry;

//      typedef boost::geometry::model::point <double, 3, cs::cartesian> point;
//      typedef boost::geometry::model::polygon <point> polygon;

//      // Get the rectange of detector points
//      double6 rect_a = a.get_image_rectangle();
//      double6 rect_b = b.get_image_rectangle();

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

    std::string type_;
    panel_list_type panel_list_;
  };

}} // namespace dxtbx::model

#endif // DXTBX_MODEL_DETECTOR_H
