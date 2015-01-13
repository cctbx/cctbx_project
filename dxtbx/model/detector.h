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
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <dxtbx/model/panel.h>
#include <dxtbx/error.h>

namespace dxtbx { namespace model {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::mat3;
  using scitbx::af::double6;
  using scitbx::af::int4;


  /**
  * A class representing a detector made up of multiple flat panel detectors.
  * The detector elements can be accessed in the same way as an array:
  *  detector[0] -> 1st detector panel.
  */
  class Detector {
  public:

    typedef std::pair<int, vec2<double> > coord_type;
    typedef Panel panel_type;
    typedef boost::ptr_vector <panel_type> panel_list_type;
    typedef panel_list_type::iterator iterator;

    /** Default constructor */
    Detector()
      : panel_list_(new panel_list_type()) {}

    /**  Construct from a single panel */
    Detector(const Panel &panel)
      : panel_list_(new panel_list_type()) {
      panel_list_->push_back(new Panel(panel));
    }

    /** Virtual destructor */
    virtual ~Detector() {}

    /** Get the begin iterator */
    iterator begin() {
      return panel_list_->begin();
    }

    /** Get the end iterator */
    iterator end() {
      return panel_list_->end();
    }

    /** Add a panel to the list of panels */
    void add_panel(Panel *panel) {
      panel_list_->push_back(panel);
    }

    /** Add a panel to the list of panels */
    panel_type& add_panel(const Panel &panel) {
      panel_list_->push_back(new Panel(panel));
      return panel_list_->back();
    }

    /** Add a panel to the list of panels */
    panel_type& add_panel() {
      return add_panel(Panel());
    }

    /** Get the number of panels */
    std::size_t size() const {
      return panel_list_->size();
    }

    /** Return a reference to a panel */
    panel_type& operator[](std::size_t index) {
      DXTBX_ASSERT(index < size());
      return (*panel_list_)[index];
    }

    /** Return a reference to a panel */
    panel_type* at(std::size_t index) {
      DXTBX_ASSERT(index < size());
      return &(*panel_list_)[index];
    }

    /** Return a const reference to a panel */
    const panel_type& operator[](std::size_t index) const {
      DXTBX_ASSERT(index < size());
      return (*panel_list_)[index];
    }

    /** Check the detector panels are the same */
    bool operator==(const Detector &detector) const {
      bool same = panel_list_->size() == detector.panel_list_->size();
      if (same) {
        for (std::size_t i = 0; i < panel_list_->size(); ++i) {
          same = same && ((*panel_list_)[i] == (*detector.panel_list_)[i]);
        }
      }
      return same;
    }

    /** Check the detector panels are not the same */
    bool operator!=(const Detector &detector) const {
      return !(*this == detector);
    }

    /** Check if the detectors are similar */
    bool is_similar_to(const Detector &rhs) const {
      bool similar = panel_list_->size() == rhs.panel_list_->size();
      if (similar) {
        for (std::size_t i = 0; i < panel_list_->size(); ++i) {
          similar = similar && ((*panel_list_)[i].is_similar_to(
            (*rhs.panel_list_)[i]));
        }
      }
      return similar;
    }

    /** Get the maximum resolution of the detector */
    double get_max_resolution(vec3<double> s0) const {
      double d_min = 0;
      for (std::size_t i = 0; i < panel_list_->size(); ++i) {
        double d = (*panel_list_)[i].get_max_resolution_at_corners(s0);
        if (d < d_min || i == 0) d_min = d;
      }
      return d_min;
    }

    /** Get ray intersection with detector */
    coord_type get_ray_intersection(vec3<double> s1) const {
      coord_type pxy(-1, vec2<double>(0, 0));
      double w_max = 0;

      // Loop through all detectors. If the w component of the (u, v, w)
      // vector points in the correct direction, and is greater than that of
      // the current closest valid coordinate, then calculate the coordinate.
      // If the coordinate is valid, then set this coordinate as the current
      // best bet.
      for (std::size_t i = 0; i < panel_list_->size(); ++i) {
        vec3 <double> v = (*panel_list_)[i].get_D_matrix() * s1;
        if (v[2] > w_max) {
          vec2<double> xy_temp(v[0] / v[2], v[1] / v[2]);
          if ((*panel_list_)[i].is_coord_valid_mm(xy_temp)) {
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

    /** finds the panel id with which s1 intersects.  Returns -1 if none do. **/
    int get_panel_intersection(vec3<double> s1) {
      int found_panel = -1;
      for (std::size_t i = 0; i < panel_list_->size(); ++i) {
        vec2<double> intersection = (*panel_list_)[i].get_ray_intersection(s1);
        if ((*panel_list_)[i].is_coord_valid_mm(intersection)) {
          DXTBX_ASSERT(found_panel < 0);
          found_panel = (int)i;
        }
      }
      return found_panel;
    }


    /** Check if any panels intersect */
//    bool do_panels_intersect() const {
//      for (std::size_t j = 0; j < panel_list_->size()-1; ++j) {
//        for (std::size_t i = j+1; i < panel_list_->size(); ++i) {
//          if (panels_intersect((*panel_list_)[j], (*panel_list_)[i])) {
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
//    panels_intersect(const Panel &a, const Panel &b) {

//      using namespace boost::geometry;

//      typedef boost::geometry::model::point <double, 3, cs::cartesian> point;
//      typedef boost::geometry::model::polygon <point> polygon;

//      // Get coordinates for panel a
//      std::size_t aw = a.get_image_size()[0];
//      std::size_t ah = a.get_image_size()[1];
//      vec3<double> a00 = a.get_pixel_lab_coord(vec2<double>(0, 0));
//      vec3<double> a01 = a.get_pixel_lab_coord(vec2<double>(aw, 0));
//      vec3<double> a10 = a.get_pixel_lab_coord(vec2<double>(0, ah));
//      vec3<double> a11 = a.get_pixel_lab_coord(vec2<double>(aw, ah));

//      // Get coordinates for panel a
//      std::size_t bw = b.get_image_size()[0];
//      std::size_t bh = b.get_image_size()[1];
//      vec3<double> b00 = b.get_pixel_lab_coord(vec2<double>(0, 0));
//      vec3<double> b01 = b.get_pixel_lab_coord(vec2<double>(bw, 0));
//      vec3<double> b10 = b.get_pixel_lab_coord(vec2<double>(0, bh));
//      vec3<double> b11 = b.get_pixel_lab_coord(vec2<double>(bw, bh));

//      // Create a polygon for the panel a plane
//      polygon poly_a;
//      append(poly_a, point(a00[0], a00[1], a00[2]));
//      append(poly_a, point(a01[0], a01[1], a01[2]));
//      append(poly_a, point(a11[0], a11[1], a11[2]));
//      append(poly_a, point(a10[0], a10[1], a10[2]));
//      append(poly_a, point(a00[0], a00[1], a00[2]));

//      // Create a polygon for the panel b plane
//      polygon poly_b;
//      append(poly_b, point(b00[0], b00[1], b00[2]));
//      append(poly_b, point(b01[0], b01[1], b01[2]));
//      append(poly_b, point(b11[0], b11[1], b11[2]));
//      append(poly_b, point(b10[0], b10[1], b10[2]));
//      append(poly_b, point(b00[0], b00[1], b00[2]));

//      // Check if the polygons intersect
//      return intersects(poly_a, poly_b);
//    }

    friend std::ostream& operator<<(std::ostream &os, const Detector &d);

    boost::shared_ptr<panel_list_type> panel_list_;
  };

  /** Print detector information */
  inline
  std::ostream& operator<<(std::ostream &os, const Detector &d) {
    os << "Detector:\n";
    for (std::size_t i = 0; i < d.size(); ++i) {
      os << d[i];
    }
    return os;
  }

}} // namespace dxtbx::model

#endif // DXTBX_MODEL_DETECTOR_H
