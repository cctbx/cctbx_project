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
#ifndef DXTBX_MODEL_DETECTOR2_H
#define DXTBX_MODEL_DETECTOR2_H

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
#include <dxtbx/model/panel2.h>
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
  class Detector2 {
  public:

    // Panel list typedefs
    typedef Panel2 panel_type;
    typedef boost::ptr_vector <panel_type> panel_list_type;
    typedef panel_list_type::iterator iterator;

    /** Default constructor */
    Detector2():
      panel_list_(new panel_list_type()) {}

    /** Virtual destructor */
    virtual ~Detector2() {}

    /** Get the begin iterator */
    iterator begin() {
      return panel_list_->begin();
    }

    /** Get the end iterator */
    iterator end() {
      return panel_list_->end();
    }

    /** Add a panel to the list of panels */
    panel_type& add_panel() {
      panel_list_->push_back(new Panel2());
      return panel_list_->back();
    }

    /** Get the number of panels */
    std::size_t size() const {
      return panel_list_->size();
    }

    /** Return a reference to a panel */
    panel_type& operator[](std::size_t index) {
      return (*panel_list_)[index];
    }

    /** Return a const reference to a panel */
    const panel_type& operator[](std::size_t index) const {
      return (*panel_list_)[index];
    }

    /** Check the detector panels are the same */
    bool operator==(const Detector2 &detector) const {
      bool same = panel_list_->size() == detector.panel_list_->size();
      if (same) {
        for (std::size_t i = 0; i < panel_list_->size(); ++i) {
          same = same && ((*panel_list_)[i] == (*detector.panel_list_)[i]);
        }
      }
      return same;
    }

    /** Check the detector panels are not the same */
    bool operator!=(const Detector2 &detector) const {
      return !(*this == detector);
    }

  protected:

    boost::shared_ptr<panel_list_type> panel_list_;
  };

}} // namespace dxtbx::model

#endif // DXTBX_MODEL_DETECTOR2_H
