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

  /**
   * A class to represent what's in an experiment.
   *
   * Contains:
   *   - imageset Access to the image data
   *   - beam The beam model
   *   - detector The detector model
   *   - goniometer The goniometer model
   *   - scan The scan model
   *   - crystal The crystal model
   *   - profile The profile model
   *   - scaling_model The scaling model
   *
   * Some of these may be set to "None"
   *
   */
  class Experiment {
  public:

    Experiment() {}

    /**
     * Initialise the experiment with models
     */
    Experiment(
          boost::shared_ptr<BeamBase> beam,
          boost::shared_ptr<Detector> detector,
          boost::shared_ptr<Goniometer> goniometer,
          boost::shared_ptr<Scan> scan,
          boost::shared_ptr<CrystalBase> crystal,
          boost::python::object profile,
          boost::python::object imageset,
          boost::python::object scaling_model)
      : beam_(beam),
        detector_(detector),
        goniometer_(goniometer),
        scan_(scan),
        crystal_(crystal),
        profile_(profile),
        imageset_(imageset),
        scaling_model_(scaling_model){}

    /**
     * Check if the beam model is the same.
     */
    bool contains(const boost::shared_ptr<BeamBase> &beam) const {
      return beam_ == beam;
    }

    /**
     * Check if the detector model is the same.
     */
    bool contains(const boost::shared_ptr<Detector> &detector) const {
      return detector_ == detector;
    }

    /**
     * Check if the detector model is the same.
     */
    bool contains(const boost::shared_ptr<Goniometer> &goniometer) const {
      return goniometer_ == goniometer;
    }

    /**
     * Check if the goniometer model is the same.
     */
    bool contains(const boost::shared_ptr<Scan> &scan) const {
      return scan_ == scan;
    }

    /**
     * Check if the crystal model is the same.
     */
    bool contains(const boost::shared_ptr<CrystalBase> &crystal) const {
      return crystal_ == crystal;
    }

    /**
     * Check models are the same.
     */
    bool contains(boost::python::object obj) const {
      boost::python::extract< boost::shared_ptr<BeamBase> > get_beam(obj);
      boost::python::extract< boost::shared_ptr<Detector> > get_detector(obj);
      boost::python::extract< boost::shared_ptr<Goniometer> > get_goniometer(obj);
      boost::python::extract< boost::shared_ptr<Scan> > get_scan(obj);
      boost::python::extract< boost::shared_ptr<CrystalBase> > get_crystal(obj);
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
      return profile_ == obj || imageset_ == obj || scaling_model_ == obj;
    }

    /**
     * Compare this experiment with another
     */
    bool operator==(const Experiment &other) const {
      return imageset_ == other.imageset_
          && beam_ == other.beam_
          && detector_ == other.detector_
          && goniometer_ == other.goniometer_
          && scan_ == other.scan_
          && profile_ == other.profile_
          && scaling_model_ == other.scaling_model_;
    }

    /**
     * Check that the experiment is consistent
     */
    bool is_consistent() const {
      return true; // FIXME
    }

    /**
     * Set the beam model
     */
    void set_beam(boost::shared_ptr<BeamBase> beam) {
      beam_ = beam;
    }

    /**
     * Get the beam model
     */
    boost::shared_ptr<BeamBase> get_beam() const {
      return beam_;
    }

    /**
     * Get the detector model
     */
    void set_detector(boost::shared_ptr<Detector> detector) {
      detector_ = detector;
    }

    /**
     * Get the detector model
     */
    boost::shared_ptr<Detector> get_detector() const {
      return detector_;
    }

    /**
     * Get the goniometer model
     */
    void set_goniometer(boost::shared_ptr<Goniometer> goniometer) {
      goniometer_ = goniometer;
    }

    /**
     * Get the goniometer model
     */
    boost::shared_ptr<Goniometer> get_goniometer() const {
      return goniometer_;
    }

    /**
     * Get the scan model
     */
    void set_scan(boost::shared_ptr<Scan> scan) {
      scan_ = scan;
    }

    /**
     * Get the scan model
     */
    boost::shared_ptr<Scan> get_scan() const {
      return scan_;
    }

    /**
     * Get the crystal model
     */
    void set_crystal(boost::shared_ptr<CrystalBase> crystal) {
      crystal_ = crystal;
    }

    /**
     * Get the crystal model
     */
    boost::shared_ptr<CrystalBase> get_crystal() const {
      return crystal_;
    }

    /**
     * Get the profile model
     */
    void set_profile(boost::python::object profile) {
      profile_ = profile;
    }

    /**
     * Get the profile model
     */
    boost::python::object get_profile() const {
      return profile_;
    }

    /**
     * Get the imageset model
     */
    void set_imageset(boost::python::object imageset) {
      imageset_ = imageset;
    }

    /**
     * Get the imageset model
     */
    boost::python::object get_imageset() const {
      return imageset_;
    }

    /**
    * Set the scaling model
    */
    void set_scaling_model(boost::python::object scaling_model) {
      scaling_model_ = scaling_model;
    }

    /**
    * Get the scaling model
    */
    boost::python::object get_scaling_model() const {
      return scaling_model_;
    }

  protected:

    boost::shared_ptr<BeamBase> beam_;
    boost::shared_ptr<Detector> detector_;
    boost::shared_ptr<Goniometer> goniometer_;
    boost::shared_ptr<Scan> scan_;
    boost::shared_ptr<CrystalBase> crystal_;
    boost::python::object profile_;
    boost::python::object imageset_;
    boost::python::object scaling_model_;

  };

}} // namespace dxtbx::model

#endif // DXTBX_MODEL_EXPERIMENT_H
