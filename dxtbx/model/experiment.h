/*
 * experiment.h
 *
 *  Copyright (C) 2017 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DXTBX_MODEL_EXPERIMENT_H
#define DXTBX_MODEL_EXPERIMENT_H

#include <iostream>
#include <cmath>
#include <boost/shared_ptr.hpp>
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <scitbx/vec3.h>
#include <scitbx/array_family/simple_io.h>
#include <scitbx/array_family/simple_tiny_io.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/scan.h>
#include <dxtbx/model/crystal.h>
#include <dxtbx/error.h>

namespace dxtbx { namespace model {

  class Experiment {
  public:

    Experiment() {}

    Experiment(
          boost::shared_ptr<Beam> beam,
          boost::shared_ptr<Detector> detector,
          boost::shared_ptr<Goniometer> goniometer,
          boost::shared_ptr<Scan> scan,
          boost::shared_ptr<Crystal> crystal,
          boost::python::object profile,
          boost::python::object imageset)
      : beam_(beam),
        detector_(detector),
        goniometer_(goniometer),
        scan_(scan),
        crystal_(crystal),
        profile_(profile),
        imageset_(imageset) {}

    bool contains(const boost::shared_ptr<Beam> &beam) const {
      return beam_ == beam;
    }

    bool contains(const boost::shared_ptr<Detector> &detector) const {
      return detector_ == detector;
    }

    bool contains(const boost::shared_ptr<Goniometer> &goniometer) const {
      return goniometer_ == goniometer;
    }

    bool contains(const boost::shared_ptr<Scan> &scan) const {
      return scan_ == scan;
    }

    bool contains(const boost::shared_ptr<Crystal> &crystal) const {
      return crystal_ == crystal;
    }

    bool contains(boost::python::object obj) const {
      boost::python::extract< boost::shared_ptr<Beam> > get_beam(obj);
      boost::python::extract< boost::shared_ptr<Detector> > get_detector(obj);
      boost::python::extract< boost::shared_ptr<Goniometer> > get_goniometer(obj);
      boost::python::extract< boost::shared_ptr<Scan> > get_scan(obj);
      boost::python::extract< boost::shared_ptr<Crystal> > get_crystal(obj);
      if (get_beam.check()) {
        return contains(get_beam());
      } else if (get_detector.check()) {
        return contains(get_detector());
      } else if (get_goniometer.check()) {
        return contains(get_goniometer());
      } else if (get_scan.check()) {
        return contains(get_scan());
      } else if (get_crystal.check()) {
        return contains(get_crystal());
      }
      return profile_ == obj || imageset_ == obj;
    }

    bool operator==(const Experiment &other) const {
      return imageset_   == other.imageset_
          && beam_       == other.beam_
          && detector_   == other.detector_
          && goniometer_ == other.goniometer_
          && scan_       == other.scan_
          && profile_    == other.profile_;
    }

    bool is_consistent() const {
      return true; // FIXME
    }

    void set_beam(boost::shared_ptr<Beam> beam) {
      beam_ = beam;
    }

    boost::shared_ptr<Beam> get_beam() const {
      return beam_;
    }

    void set_detector(boost::shared_ptr<Detector> detector) {
      detector_ = detector;
    }

    boost::shared_ptr<Detector> get_detector() const {
      return detector_;
    }

    void set_goniometer(boost::shared_ptr<Goniometer> goniometer) {
      goniometer_ = goniometer;
    }

    boost::shared_ptr<Goniometer> get_goniometer() const {
      return goniometer_;
    }

    void set_scan(boost::shared_ptr<Scan> scan) {
      scan_ = scan;
    }

    boost::shared_ptr<Scan> get_scan() const {
      return scan_;
    }

    void set_crystal(boost::shared_ptr<Crystal> crystal) {
      crystal_ = crystal;
    }

    boost::shared_ptr<Crystal> get_crystal() const {
      return crystal_;
    }

    void set_profile(boost::python::object profile) {
      profile_ = profile;
    }

    boost::python::object get_profile() const {
      return profile_;
    }

    void set_imageset(boost::python::object imageset) {
      imageset_ = imageset;
    }

    boost::python::object get_imageset() const {
      return imageset_;
    }

  protected:

    boost::shared_ptr<Beam> beam_;
    boost::shared_ptr<Detector> detector_;
    boost::shared_ptr<Goniometer> goniometer_;
    boost::shared_ptr<Scan> scan_;
    boost::shared_ptr<Crystal> crystal_;
    boost::python::object profile_;
    boost::python::object imageset_;

  };

}} // namespace dxtbx::model

#endif // DXTBX_MODEL_EXPERIMENT_H
