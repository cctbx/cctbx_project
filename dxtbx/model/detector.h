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

  /**
  * A class representing a detector made up of multiple flat panel detectors.
  * The detector elements can be accessed in the same way as an array:
  *  detector[0] -> 1st detector panel.
  */
  class Detector {
  public:

    typedef std::pair<int, vec2<double> > coord_type;

    // Panel list typedefs
    typedef Panel panel_type;
    typedef scitbx::af::shared <panel_type> panel_list_type;
    typedef panel_list_type::iterator iterator;

    /** Default constructor */
    Detector() {}

    /** Initialise the detector with singla panel */
    Detector(const Panel &panel) {
      panel_list_.push_back(panel);
    }

    /** Initialise the detector with an array of panels */
    Detector(const panel_list_type &panel_list)
      : panel_list_(panel_list) {}

    /** Virtual destructor */
    virtual ~Detector() {}

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
    bool operator==(const Detector &detector) {
      bool same = panel_list_.size() == detector.panel_list_.size();
      if (same) {
        for (std::size_t i = 0; i < panel_list_.size(); ++i) {
          same = same && (panel_list_[i] == detector.panel_list_[i]);
        }
      }
      return same;
    }

    /** Check the detector panels are not the same */
    bool operator!=(const Detector &detector) {
      return !(*this == detector);
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
    double get_max_resolution(vec3<double> s0, double wavelength) {
      double d_min = 0;
      for (std::size_t i = 0; i < panel_list_.size(); ++i) {
        double d = panel_list_[i].get_max_resolution_at_corners(s0, wavelength);
        if (d < d_min || i == 0) d_min = d;
      }
      return d_min;
    }

    /** Get ray intersection with detector */
    coord_type get_ray_intersection(vec3<double> s1) {
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
          if (panel_list_[i].is_coord_valid(xy_temp) && v[2] > w_max) {
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

    friend std::ostream& operator<< (std::ostream& , const Detector&);

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

    panel_list_type panel_list_;
  };

  /** Print the detector information to the ostream */
  inline
  std::ostream& operator<< (std::ostream &os, const Detector &d) {
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

}} // namespace dxtbx::model

#endif // DXTBX_MODEL_DETECTOR_H
