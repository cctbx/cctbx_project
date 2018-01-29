/*
 * experiment_list.h
 *
 *  Copyright (C) 2017 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DXTBX_MODEL_EXPERIMENT_LIST_H
#define DXTBX_MODEL_EXPERIMENT_LIST_H

#include <iostream>
#include <cmath>
#include <boost/shared_ptr.hpp>
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <scitbx/vec3.h>
#include <scitbx/array_family/simple_io.h>
#include <scitbx/array_family/simple_tiny_io.h>
#include <dxtbx/model/experiment.h>
#include <dxtbx/error.h>

namespace dxtbx { namespace model {

  /**
   * This class contains a list of experiments
   */
  class ExperimentList {
  public:

    typedef scitbx::af::shared<Experiment> shared_type;
    typedef scitbx::af::const_ref<Experiment> const_ref_type;
    typedef shared_type::const_iterator const_iterator;
    typedef shared_type::iterator iterator;

    ExperimentList() {}

    /**
     * Initialize from the data
     */
    ExperimentList(const const_ref_type &data)
      : data_(data.begin(), data.end()) {}

    /**
     * Get a shared array of experiments
     */
    shared_type data() const {
      return data_;
    }

    /**
     * Get the experiment at the given index
     */
    Experiment& operator[](std::size_t index) {
      DXTBX_ASSERT(index < data_.size());
      return data_[index];
    }

    /**
     * Get the experiment at the given index
     */
    const Experiment& operator[](std::size_t index) const {
      DXTBX_ASSERT(index < data_.size());
      return data_[index];
    }

    /**
     * Erase the experiment at the given index
     */
    void erase(std::size_t index) {
      DXTBX_ASSERT(index < data_.size());
      data_.erase(data_.begin()+index, data_.begin()+index+1);
    }

    /**
     * Get the number of experiments
     */
    std::size_t size() const {
      return data_.size();
    }

    /**
     * Get an iterator to the beginning of the list
     */
    const_iterator begin() const {
      return data_.begin();
    }

    /**
     * Get an iterator to the end of the list
     */
    const_iterator end() const {
      return data_.end();
    }

    /**
     * Get an iterator to the beginning of the list
     */
    iterator begin() {
      return data_.begin();
    }

    /**
     * Get an iterator to the end of the list
     */
    iterator end() {
      return data_.end();
    }

    /**
     * Clear the experiments
     */
    void clear() {
      data_.clear();
    }

    /**
     * Check if the experiment list is empty
     */
    bool empty() const {
      return data_.empty();
    }

    /**
     * Append an experiment to the list
     */
    void append(const Experiment &experiment) {
      data_.push_back(experiment);
    }

    /**
     * Extend the experiment list with experiments from another
     */
    void extend(const ExperimentList &experiment_list) {
      data_.insert(
          data_.end(),
          experiment_list.data_.begin(),
          experiment_list.data_.end());
    }

    /**
     * Check if an experiment contains the beam model
     */
    bool contains(const boost::shared_ptr<BeamBase> &beam) const {
      for (std::size_t i = 0; i < size(); ++i) {
        if (data_[i].contains(beam)) {
          return true;
        }
      }
      return false;
    }

    /**
     * Check if an experiment contains the detector model
     */
    bool contains(const boost::shared_ptr<Detector> &detector) const {
      for (std::size_t i = 0; i < size(); ++i) {
        if (data_[i].contains(detector)) {
          return true;
        }
      }
      return false;
    }

    /**
     * Check if an experiment contains the goniometer model
     */
    bool contains(const boost::shared_ptr<Goniometer> &goniometer) const {
      for (std::size_t i = 0; i < size(); ++i) {
        if (data_[i].contains(goniometer)) {
          return true;
        }
      }
      return false;
    }

    /**
     * Check if an experiment contains the scan model
     */
    bool contains(const boost::shared_ptr<Scan> &scan) const {
      for (std::size_t i = 0; i < size(); ++i) {
        if (data_[i].contains(scan)) {
          return true;
        }
      }
      return false;
    }

    /**
     * Check if an experiment contains the crystal model
     */
    bool contains(const boost::shared_ptr<CrystalBase> &crystal) const {
      for (std::size_t i = 0; i < size(); ++i) {
        if (data_[i].contains(crystal)) {
          return true;
        }
      }
      return false;
    }

    /**
     * Check if an experiment contains the model
     */
    bool contains(boost::python::object obj) const {
      for (std::size_t i = 0; i < size(); ++i) {
        if (data_[i].contains(obj)) {
          return true;
        }
      }
      return false;
    }

    /**
     * Replace all beam models
     */
    void replace(
        boost::shared_ptr<BeamBase> a,
        boost::shared_ptr<BeamBase> b) {
      for (std::size_t i = 0; i < size(); ++i) {
        if (data_[i].get_beam() == a) {
          data_[i].set_beam(b);
        }
      }
    }

    /**
     * Replace all detector models
     */
    void replace(
        boost::shared_ptr<Detector> a,
        boost::shared_ptr<Detector> b) {
      for (std::size_t i = 0; i < size(); ++i) {
        if (data_[i].get_detector() == a) {
          data_[i].set_detector(b);
        }
      }
    }

    /**
     * Replace all goniometer models
     */
    void replace(
        boost::shared_ptr<Goniometer> a,
        boost::shared_ptr<Goniometer> b) {
      for (std::size_t i = 0; i < size(); ++i) {
        if (data_[i].get_goniometer() == a) {
          data_[i].set_goniometer(b);
        }
      }
    }

    /**
     * Replace all scan models
     */
    void replace(
        boost::shared_ptr<Scan> a,
        boost::shared_ptr<Scan> b) {
      for (std::size_t i = 0; i < size(); ++i) {
        if (data_[i].get_scan() == a) {
          data_[i].set_scan(b);
        }
      }
    }

    /**
     * Replace all crystal models
     */
    void replace(
        boost::shared_ptr<CrystalBase> a,
        boost::shared_ptr<CrystalBase> b) {
      for (std::size_t i = 0; i < size(); ++i) {
        if (data_[i].get_crystal() == a) {
          data_[i].set_crystal(b);
        }
      }
    }

    /**
     * Replace all other models
     */
    void replace(
        boost::python::object a,
        boost::python::object b) {
      boost::python::extract< boost::shared_ptr<BeamBase> > get_beam_a(a);
      boost::python::extract< boost::shared_ptr<BeamBase> > get_beam_b(b);
      boost::python::extract< boost::shared_ptr<Detector> > get_detector_a(a);
      boost::python::extract< boost::shared_ptr<Detector> > get_detector_b(b);
      boost::python::extract< boost::shared_ptr<Goniometer> > get_goniometer_a(a);
      boost::python::extract< boost::shared_ptr<Goniometer> > get_goniometer_b(b);
      boost::python::extract< boost::shared_ptr<Scan> > get_scan_a(a);
      boost::python::extract< boost::shared_ptr<Scan> > get_scan_b(b);
      boost::python::extract< boost::shared_ptr<CrystalBase> > get_crystal_a(a);
      boost::python::extract< boost::shared_ptr<CrystalBase> > get_crystal_b(b);
      if (get_beam_a.check()) {
        DXTBX_ASSERT(get_beam_b.check());
        replace(get_beam_a(), get_beam_b());
      } else if (get_detector_a.check()) {
        DXTBX_ASSERT(get_detector_b.check());
        replace(get_detector_a(), get_detector_b());
      } else if (get_goniometer_a.check()) {
        DXTBX_ASSERT(get_goniometer_b.check());
        replace(get_goniometer_a(), get_goniometer_b());
      } else if (get_scan_a.check()) {
        DXTBX_ASSERT(get_scan_b.check());
        replace(get_scan_a(), get_scan_b());
      } else if (get_crystal_a.check()) {
        DXTBX_ASSERT(get_crystal_b.check());
        replace(get_crystal_a(), get_crystal_b());
      } else {
        DXTBX_ASSERT(!get_beam_b.check());
        DXTBX_ASSERT(!get_detector_b.check());
        DXTBX_ASSERT(!get_goniometer_b.check());
        DXTBX_ASSERT(!get_scan_b.check());
        DXTBX_ASSERT(!get_crystal_b.check());
        for (std::size_t i = 0; i < size(); ++i) {
          if (data_[i].get_profile() == a) {
            data_[i].set_profile(b);
          }
          if (data_[i].get_imageset() == a) {
            data_[i].set_imageset(b);
          }
        }
      }
    }

    /**
     * Get indices which have this model
     */
    scitbx::af::shared<std::size_t> indices(const boost::shared_ptr<BeamBase> &obj) const {
      scitbx::af::shared<std::size_t> result;
      for (std::size_t i = 0; i < size(); ++i) {
         if (data_[i].get_beam() == obj) {
             result.push_back(i);
         }
      }
      return result;
    }

    /**
     * Get indices which have this model
     */
    scitbx::af::shared<std::size_t> indices(const boost::shared_ptr<Detector> &obj) const {
      scitbx::af::shared<std::size_t> result;
      for (std::size_t i = 0; i < size(); ++i) {
         if (data_[i].get_detector() == obj) {
             result.push_back(i);
         }
      }
      return result;
    }

    /**
     * Get indices which have this model
     */
    scitbx::af::shared<std::size_t> indices(const boost::shared_ptr<Goniometer> &obj) const {
      scitbx::af::shared<std::size_t> result;
      for (std::size_t i = 0; i < size(); ++i) {
         if (data_[i].get_goniometer() == obj) {
             result.push_back(i);
         }
      }
      return result;
    }

    /**
     * Get indices which have this model
     */
    scitbx::af::shared<std::size_t> indices(const boost::shared_ptr<Scan> &obj) const {
      scitbx::af::shared<std::size_t> result;
      for (std::size_t i = 0; i < size(); ++i) {
         if (data_[i].get_scan() == obj) {
             result.push_back(i);
         }
      }
      return result;
    }

    /**
     * Get indices which have this model
     */
    scitbx::af::shared<std::size_t> indices(const boost::shared_ptr<CrystalBase> &obj) const {
      scitbx::af::shared<std::size_t> result;
      for (std::size_t i = 0; i < size(); ++i) {
         if (data_[i].get_crystal() == obj) {
             result.push_back(i);
         }
      }
      return result;
    }

    /**
     * Get indices which have this model
     */
    scitbx::af::shared<std::size_t> indices(boost::python::object obj) const {
      boost::python::extract< boost::shared_ptr<BeamBase> > get_beam(obj);
      boost::python::extract< boost::shared_ptr<Detector> > get_detector(obj);
      boost::python::extract< boost::shared_ptr<Goniometer> > get_goniometer(obj);
      boost::python::extract< boost::shared_ptr<Scan> > get_scan(obj);
      boost::python::extract< boost::shared_ptr<CrystalBase> > get_crystal(obj);
      if (get_beam.check()) {
        return indices(get_beam());
      } else if (get_detector.check()) {
        return indices(get_detector());
      } else if (get_goniometer.check()) {
        return indices(get_goniometer());
      } else if (get_scan.check()) {
        return indices(get_scan());
      } else if (get_crystal.check()) {
        return indices(get_crystal());
      }
      scitbx::af::shared<std::size_t> result;
      for (std::size_t i = 0; i < size(); ++i) {
        if (data_[i].get_profile() == obj ||
            data_[i].get_imageset() == obj ||
            data_[i].get_scaling_model() == obj) {
          result.push_back(i);
        }
      }
      return result;
    }

    /**
     * Get indices which match all the models
     */
    scitbx::af::shared<std::size_t> where(
        boost::shared_ptr<BeamBase> beam,
        boost::shared_ptr<Detector> detector,
        boost::shared_ptr<Goniometer> goniometer,
        boost::shared_ptr<Scan> scan,
        boost::shared_ptr<CrystalBase> crystal,
        boost::python::object profile,
        boost::python::object imageset,
        boost::python::object scaling_model) const {
      scitbx::af::shared<std::size_t> result;
      for (std::size_t i = 0; i < size(); ++i) {
        if (beam && data_[i].get_beam() != beam) {
          continue;
        }
        if (detector && data_[i].get_detector() != detector) {
          continue;
        }
        if (goniometer && data_[i].get_goniometer() != goniometer) {
          continue;
        }
        if (scan && data_[i].get_scan() != scan) {
          continue;
        }
        if (crystal && data_[i].get_crystal() != crystal) {
          continue;
        }
        if (profile && data_[i].get_profile() != profile) {
          continue;
        }
        if (imageset && data_[i].get_imageset() != imageset) {
          continue;
        }
        if (scaling_model && data_[i].get_scaling_model() != scaling_model){
          continue;
        }
        result.push_back(i);
      }
      return result;
    }

    /**
     * Check if experiments are consistent
     */
    bool is_consistent() const {
      for (std::size_t i = 0; i < size(); ++i) {
        if (!data_[i].is_consistent()) {
          return false;
        }
      }
      return true;
    }

  protected:

    shared_type data_;
  };

}} // namespace dxtbx::model

#endif // DXTBX_MODEL_EXPERIMENT_LIST_H
