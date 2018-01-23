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
#include <vector>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/iterator/indirect_iterator.hpp>
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


  namespace detail {

    /**
     * @return True/False is point p of the left of the line between a and b
     */
    inline
    bool is_left(const vec2<double> &a, const vec2<double> &b, const vec2<double> &p) {
      return ((b[0] - a[0])*(p[1] - a[1]) - (b[1] - a[1])*(p[0] - a[0])) > 0;
    }

    /**
     * Compute the convex hull using the "gift wrapping" algorithm
     * @param x The set of input points
     * @return The points in the convex hull
     */
    inline
    scitbx::af::shared< vec2<double> > convex_hull(const scitbx::af::const_ref< vec2<double> > &x) {

      DXTBX_ASSERT(x.size() > 2);

      scitbx::af::shared< vec2<double> > result;

      // Find the leftmost point
      std::size_t current = 0;
      for (std::size_t i = 1 ; i < x.size(); ++i) {
        if (x[i][0] < x[current][0]) {
          current = i;
        }
      }

      // Save the starting index
      std::size_t first = current;
      while (true) {

        // Add the current point to the results
        result.push_back(x[current]);

        // Find the next point in the convex hull
        std::size_t last = 0;
        for (std::size_t i = 1; i < x.size(); ++i) {
          if (last == current || is_left(x[current], x[last], x[i])) {
            last = i;
          }
        }
        current = last;

        // Break the loop if we're back at the start
        if (last == first) {
          break;
        }
      }

      // Return the convex hull
      return result;
    }

    /**
     * Compute the distance from a point P to a line segment defined by points a
     * and b
     */
    inline
    double distance_to_line_segment(const vec2<double> &a, const vec2<double> &b, const vec2<double> &p) {
      vec2<double> l = (b - a);
      double t = l * (p - a) / l.length_sq();
      double distance = 0.0;
      if (t <= 0) {
        distance = (p - a).length();
      } else if (t >= 1.0) {
        distance = (p - (a + l)).length();
      } else {
        distance = (p - (a + t * l)).length();
      }
      return distance;
    }

  }

  /**
   * A class representing a multi-panel hierarchical detector.
   *
   * The detector is a composite object made up of 1 or more panels and optionally
   * a hierarchy which encodes groups of panels. The panels can be accessed as a
   * flat data structure (i.e. an array) and the panels and groups can be accessed
   * through a tree interface.
   */
  class Detector {
  public:



    /**
     * A class representing a node in the detector hierarchy. The node class
     * inherits from the panel class and adds some methods for adding children,
     * querying parents and roots etc.
     *
     */
    class Node : public Panel {
    public:

      typedef Node* pointer;
      typedef const Node* const_pointer;
      typedef boost::ptr_vector<Node>::iterator iterator;
      typedef boost::ptr_vector<Node>::const_iterator const_iterator;

      /**
       * Construct using a parent detector reference. The detector reference keeps
       * hold of the flat array of panels so that when the hierarchy is built, the
       * panels are added automatically to the flat array of panels.
       */
      Node(Detector *detector)
        : detector_(detector),
          parent_(NULL),
          is_panel_(false) {}

      /**
       * Construct using a parent detector reference. The detector reference keeps
       * hold of the flat array of panels so that when the hierarchy is built, the
       * panels are added automatically to the flat array of panels.
       */
      Node(Detector *detector, const Panel &panel)
        : Panel(panel),
          detector_(detector),
          parent_(NULL),
          is_panel_(false) {}

      /**
       * Add a group to the detector node.
       */
      pointer add_group() {
        DXTBX_ASSERT(!is_panel());
        Node *node = new Node(detector_);
        node->parent_ = this;
        node->is_panel_ = false;
        node->set_parent_frame(
            get_fast_axis(),
            get_slow_axis(),
            get_origin());
        children_.push_back(node);
        return node;
      }

      /**
       * Add a group to the detector node.
       */
      pointer add_group(const Panel &group) {
        DXTBX_ASSERT(!is_panel());
        Node *node = new Node(detector_, group);
        node->parent_ = this;
        node->is_panel_ = false;
        node->set_parent_frame(
            get_fast_axis(),
            get_slow_axis(),
            get_origin());
        children_.push_back(node);
        return node;
      }

      /**
       * Add a panel to the detector node
       */
      pointer add_panel() {
        DXTBX_ASSERT(!is_panel());
        Node *node = new Node(detector_);
        node->parent_ = this;
        node->is_panel_ = true;
        node->set_parent_frame(
            get_fast_axis(),
            get_slow_axis(),
            get_origin());
        children_.push_back(node);
        detector_->data_->panels.push_back(node);
        return node;
      }

      /**
       * Add a panel to the detector node
       */
      pointer add_panel(const Panel &panel) {
        DXTBX_ASSERT(!is_panel());
        Node *node = new Node(detector_, panel);
        node->parent_ = this;
        node->is_panel_ = true;
        node->set_parent_frame(
            get_fast_axis(),
            get_slow_axis(),
            get_origin());
        children_.push_back(node);
        detector_->data_->panels.push_back(node);
        return node;
      }

      /**
       * Add a panel to the detector node
       */
      pointer add_panel(const Panel &panel, std::size_t index) {
        DXTBX_ASSERT(!is_panel());
        Node *node = new Node(detector_, panel);
        node->parent_ = this;
        node->is_panel_ = true;
        node->set_parent_frame(
            get_fast_axis(),
            get_slow_axis(),
            get_origin());
        children_.push_back(node);
        if (detector_->data_->panels.size() <= index) {
          detector_->data_->panels.resize(index+1, NULL);
        }
        DXTBX_ASSERT(detector_->data_->panels[index] == NULL);
        detector_->data_->panels[index] = node;
        return node;
      }

      /**
       * Get the node's parent
       */
      const_pointer parent() const {
        return parent_;
      }

      /**
       * Get the node's parent
       */
      pointer parent() {
        return parent_;
      }

      /**
       * Get the tree root
       */
      const_pointer root() const {
        return parent_ == NULL ? this : parent_->root();
      }

      /**
       * Get the tree root
       */
      pointer root() {
        return parent_ == NULL ? this : parent_->root();
      }

      /**
       * Get the child node at the given index
       */
      const_pointer operator[](std::size_t i) const {
        DXTBX_ASSERT(i < children_.size());
        return &children_[i];
      }

      /**
       * Get the child node at the given index
       */
      pointer operator[](std::size_t i) {
        DXTBX_ASSERT(i < children_.size());
        return &children_[i];
      }

      /**
       * Get the iterator to the beginning of the child list
       */
      const_iterator begin() const {
        return children_.begin();
      }

      /**
       * Get the iterator to the beginning of the child list
       */
      iterator begin() {
        return children_.begin();
      }

      /**
       * Get the iterator to the end of the child list
       */
      const_iterator end() const {
        return children_.end();
      }

      /**
       * Get the iterator to the end of the child list
       */
      iterator end() {
        return children_.end();
      }

      /**
       * Check if the child list is empty
       */
      bool empty() const {
        return children_.empty();
      }

      /**
       * Get the number of children
       */
      std::size_t size() const {
        return children_.size();
      }

      std::size_t index() const {
        DXTBX_ASSERT(is_panel());
        for (std::size_t i = 0; i < detector_->size(); ++i) {
          if (this == detector_->data_->panels[i]) {
            return i;
          }
        }
        throw DXTBX_ERROR("Programmer Error: no panel in detector");
        return 0;
      }

      /**
       * Is the node a panel
       */
      bool is_panel() const {
        return is_panel_;
      }

      /**
       * Is the node a group
       */
      bool is_group() const {
        return !is_panel();
      }

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
        Panel::set_frame(d1, d2, d0);
        for (std::size_t i = 0; i < children_.size(); ++i) {
          children_[i].set_parent_frame(
              get_fast_axis(),
              get_slow_axis(),
              get_origin());
        }
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
        Panel::set_local_frame(d1, d2, d0);
        for (std::size_t i = 0; i < children_.size(); ++i) {
          children_[i].set_parent_frame(
              get_fast_axis(),
              get_slow_axis(),
              get_origin());
        }
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
        Panel::set_parent_frame(d1, d2, d0);
        for (std::size_t i = 0; i < children_.size(); ++i) {
          children_[i].set_parent_frame(
              get_fast_axis(),
              get_slow_axis(),
              get_origin());
        }
      }

      /**
       * Test if everything is equal
       */
      bool operator==(const Node &rhs) const {
        if (Panel::operator==(rhs)) {
          if (size() == rhs.size()) {
            bool all_eq = true;
            for (std::size_t i = 0; i < size(); ++i) {
              if (children_[i] != rhs.children_[i]) {
                all_eq = false;
                break;
              }
            }
            return all_eq;
          }
        }
        return false;
      }

      /**
       * Test inequality
       */
      bool operator!=(const Node &rhs) const {
        return !(*this == rhs);
      }

      /**
       * Check that the nodes are similar
       */
      bool is_similar_to(
          const Node &rhs,
          double fast_axis_tolerance,
          double slow_axis_tolerance,
          double origin_tolerance,
          bool static_only,
          bool ignore_trusted_range=false) const {
        if (Panel::is_similar_to(
              rhs,
              fast_axis_tolerance,
              slow_axis_tolerance,
              origin_tolerance,
              static_only,
              ignore_trusted_range)) {
          if (size() == rhs.size()) {
            bool all_eq = true;
            for (std::size_t i = 0; i < size(); ++i) {
              if (!children_[i].is_similar_to(
                    rhs.children_[i],
                    fast_axis_tolerance,
                    slow_axis_tolerance,
                    origin_tolerance,
                    static_only,
                    ignore_trusted_range)) {
                all_eq = false;
                break;
              }
            }
            return all_eq;
          }
        }
        return false;
      }

    protected:

      Detector *detector_;
      pointer parent_;
      boost::ptr_vector<Node> children_;
      bool is_panel_;
    };



    typedef std::pair<int, vec2<double> > coord_type;
    typedef Node::pointer node_pointer;
    typedef Node::const_pointer const_node_pointer;
    typedef Panel panel_type;
    typedef boost::indirect_iterator<std::vector<Node::pointer>::iterator> iterator;
    typedef boost::indirect_iterator<std::vector<Node::pointer>::const_iterator> const_iterator;

    /**
     * A helper class to contain the data
     */
    class DetectorData {
    public:

      DetectorData(Detector *detector)
        : root(detector) {}

      DetectorData(Detector *detector, const Panel &panel)
        : root(detector, panel) {}

      Node root;
      std::vector<Node::pointer> panels;
    };

    /**
     * Initialise the detector
     */
    Detector()
      : data_(boost::make_shared<DetectorData>(this)) {}

    /**
     * Copy another detector
     */
    Detector(const Detector &other)
        : data_(boost::make_shared<DetectorData>(this, *(other.root()))) {
      // The initializer copies the main panel data; now do the rest
      copy_node_subtree(root(), other.root());
      // Validate that everything appears to have been copied
      DXTBX_ASSERT(size() == other.size());
      for (std::size_t i = 0; i < size(); ++i) {
        DXTBX_ASSERT(at(i) != NULL);
      }
    }

    /**
     * Construct with a single panel
     */
    Detector(const Panel &panel)
      : data_(boost::make_shared<DetectorData>(this)) {
      add_panel(panel);
    }

    /**
     * Add a group to the root node.
     */
    node_pointer add_group() {
      return data_->root.add_group();
    }

    /**
     * Add a group to the root node.
     */
    node_pointer add_group(const Panel &group) {
      return data_->root.add_group(group);
    }

    /**
     * Add a panel to the root node
     */
    node_pointer add_panel() {
      return data_->root.add_panel();
    }

    /**
     * Add a panel to the root node
     */
    node_pointer add_panel(const Panel &panel) {
      return data_->root.add_panel(panel);
    }

    /**
     * Get the root node
     */
    node_pointer root() {
      return &(data_->root);
    }

    /**
     * Get the root node
     */
    const_node_pointer root() const {
      return &(data_->root);
    }

    /**
     * Get the panel at the current position
     */
    panel_type& operator[](std::size_t i) {
      DXTBX_ASSERT(i < data_->panels.size());
      return *(data_->panels[i]);
    }

    /** Return a const reference to a panel */
    const panel_type& operator[](std::size_t i) const {
     DXTBX_ASSERT(i < size());
      return *(data_->panels[i]);
    }

    /** Return a pointer to a panel */
    panel_type* at(std::size_t i) {
      DXTBX_ASSERT(i < size());
      return data_->panels[i];
    }

    std::size_t size() const {
      return data_->panels.size();
    }

    /** Get the begin iterator */
    const_iterator begin() const {
      return const_iterator(data_->panels.begin());
    }

    /** Get the begin iterator */
    iterator begin() {
      return iterator(data_->panels.begin());
    }

    /** Get the end iterator */
    const_iterator end() const {
      return const_iterator(data_->panels.end());
    }

    /** Get the end iterator */
    iterator end() {
      return iterator(data_->panels.end());
    }

    /** Check the detector panels are the same */
    bool operator==(const Detector &detector) const {
      bool same = false;
      if (size() == detector.size() && (data_->root == detector.data_->root)) {
        same = true;
        for (std::size_t i = 0; i < size(); ++i) {
          same = same && *(data_->panels[i]) == *(detector.data_->panels[i]);
        }
      }
      return same;
    }

    /** Check the detector panels are not the same */
    bool operator!=(const Detector &detector) const {
      return !(*this == detector);
    }

    /** Check if the detectors are similar */
    bool is_similar_to(const Detector &rhs,
                        double fast_axis_tolerance,
                        double slow_axis_tolerance,
                        double origin_tolerance,
                        bool static_only,
                        bool ignore_trusted_range=false) const {
      bool similar = false;
      if (data_->root.is_similar_to(
            rhs.data_->root,
            fast_axis_tolerance,
            slow_axis_tolerance,
            origin_tolerance,
            static_only,
            ignore_trusted_range) && size() == rhs.size()) {
        similar = true;
        for (std::size_t i = 0; i < size(); ++i) {
          similar = similar && (data_->panels[i]->is_similar_to(
            *rhs.data_->panels[i],
            fast_axis_tolerance,
            slow_axis_tolerance,
            origin_tolerance,
            static_only,
            ignore_trusted_range));
        }
      }
      return similar;
    }

    /** Get the maximum resolution of the detector */
    double get_max_resolution(vec3<double> s0) const {
      double d_min = 0;
      for (std::size_t i = 0; i < size(); ++i) {
        double d = (*this)[i].get_max_resolution_at_corners(s0);
        if (d < d_min || i == 0) d_min = d;
      }
      return d_min;
    }

    /**
     * Compute a sensible "inscribed" resolution that works in the case of
     * multiple panel detectors. This function does a stereographic projection
     * of the corners of each panel in the detector and then computes the
     * convex hull around the resulting points. Then for each line segment in
     * the convex hull, the shortest distance to zero is computed. The smallest
     * distance is then taken and converted back to an angle to compute the
     * maximum resolution.
     * @param s0 The incident beam vector
     * @returns The maximum resolution
     */
    double get_max_inscribed_resolution(vec3<double> s0) const {

      // If only 1 panel then use simple method
      try {
        if (size() == 1) {
          return (*this)[0].get_max_resolution_ellipse(s0);
        }
      } catch (dxtbx::error) {
        // do nothing
      }

      // Save the length of the s0 vector
      double s0_length = s0.length();

      // Choose some axes orthogonal s0
      vec3<double> za = -s0.normalize();
      vec3<double> ya = ((za * vec3<double>(1,0,0) < 0.9)
        ? za.cross(vec3<double>(1,0,0))
        : za.cross(vec3<double>(0,1,0))).normalize();
      vec3<double> xa = za.cross(ya).normalize();

      // Compute the stereographic projection of panel corners
      scitbx::af::shared< vec2<double> > points;
      for (std::size_t i = 0; i < size(); ++i) {
        std::size_t width = (*this)[i].get_image_size()[0];
        std::size_t height = (*this)[i].get_image_size()[1];
        scitbx::af::tiny< vec2<double>, 4 > corners;
        corners[0] = vec2<double>(0,0);
        corners[1] = vec2<double>(width,0);
        corners[2] = vec2<double>(0,height);
        corners[3] = vec2<double>(width,height);
        for (std::size_t j = 0; j < 4; ++j) {
          vec3<double> p = (*this)[i].get_pixel_lab_coord(corners[j]).normalize();
          double pdotx = p * xa;
          double pdoty = p * ya;
          double pdotz = p * za;
          points.push_back(
              vec2<double>(
                2.0 * pdotx / (1.0 - pdotz),
                2.0 * pdoty / (1.0 - pdotz)));
        }
      }

      // Compute the convex hull of points
      scitbx::af::shared< vec2<double> > hull = detail::convex_hull(points.const_ref());
      DXTBX_ASSERT(hull.size() >= 4);

      // Compute the minimum distance to the line segments
      double min_distance = -1;
      std::size_t first = hull.size() - 1;
      for (std::size_t last = 0; last < hull.size(); ++last) {
        vec2<double> zero(0,0);
        double d = detail::distance_to_line_segment(hull[first], hull[last], zero);
        if (min_distance < 0 || d < min_distance) {
          min_distance = d;
        }
        first = last;
      }
      DXTBX_ASSERT(min_distance > 0);

      // Compute and return the resolution
      double angle = 2.0 * std::atan(min_distance / 2.0);
      double den = (2.0 * s0_length * std::sin(0.5 * angle));
      DXTBX_ASSERT(den != 0);
      return 1.0 / den;
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
      for (std::size_t i = 0; i < size(); ++i) {
        vec3 <double> v = (*this)[i].get_D_matrix() * s1;
        if (v[2] > w_max) {
          vec2<double> xy_temp(v[0] / v[2], v[1] / v[2]);
          if ((*this)[i].is_coord_valid_mm(xy_temp)) {
            pxy = coord_type((int)i, xy_temp);
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
      for (std::size_t i = 0; i < size(); ++i) {
        try {
          vec2<double> intersection = (*this)[i].get_ray_intersection(s1);
          if ((*this)[i].is_coord_valid_mm(intersection)) {
            found_panel = (int)i;
            break;
          }
        } catch(dxtbx::error) {
          // pass
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

    /** Rotate the detector about an axis */
    void rotate_around_origin(vec3<double> axis, double angle) {
      for (std::size_t i = 0; i < size(); ++i) {
        (*this)[i].rotate_around_origin(axis, angle);
      }
    }


  protected:
    friend std::ostream& operator<<(std::ostream &os, const Detector &d);

    boost::shared_ptr<DetectorData> data_;

  private:
   /**
    * Copy the child panels and groups, recursively, of a detector node.
    *
    * Note that this doesn't touch the Panel/detector contents of the
    * destination node.
    *
    * Panel children are assumed to be leafs without further hierarchy.
    */
    void copy_node_subtree(Node::pointer dest, Node::const_pointer source) {
     for (Node::const_iterator it = source->begin(); it != source->end(); ++it) {
       if (it->is_panel()) {
         dest->add_panel(*it, it->index());
       } else {
         copy_node_subtree(dest->add_group(*it), &*it);
       }
     }
    }
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
