/*
 * beam.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DXTBX_MODEL_BEAM_H
#define DXTBX_MODEL_BEAM_H

#include <iostream>
#include <cmath>
#include <scitbx/vec3.h>
#include <scitbx/array_family/simple_io.h>
#include <scitbx/array_family/simple_tiny_io.h>
#include <dxtbx/error.h>
#include "model_helpers.h"

namespace dxtbx { namespace model {

  using scitbx::vec3;

  /** Base class for beam objects */
  class BeamBase {};

  /** A class to represent a simple beam. */
  class Beam : public BeamBase {
  public:
    /** Default constructor: initialise all to zero */
    Beam()
      : wavelength_(0.0),
        direction_(0.0, 0.0, 0.0) {}

    /**
     * Initialise all the beam parameters.
     * @param direction The beam direction vector.
     */
    Beam(vec3 <double> s0) {
      DXTBX_ASSERT(s0.length() > 0);
      wavelength_ = 1.0 / s0.length();
      direction_ = s0.normalize();    
    }

    /**
     * Initialise all the beam parameters. Normalize the direction vector
     * and give it the length of 1.0 / wavelength
     * @param wavelength The wavelength of the beam
     * @param direction The beam direction vector.
     */
    Beam(vec3 <double> direction, double wavelength)
      : wavelength_(wavelength) {      
      DXTBX_ASSERT(direction.length() > 0);
      direction_ = direction.normalize();    
    }

    /** Virtual destructor */
    virtual ~Beam() {}

    /** Get the direction */
    vec3 <double> get_direction() const {
      return direction_;
    }

    /** Get the wavelength */
    double get_wavelength() const {
      return wavelength_;
    }

    /** Set the direction. */
    void set_direction(vec3 <double> direction) {
      DXTBX_ASSERT(direction.length() > 0);
      direction_ = direction.normalize();
    }

    /** Set the wavelength */
    void set_wavelength(double wavelength) {
      wavelength_ = wavelength;
    }

    /** Get the wave vector in units of inverse angstroms */
    vec3 <double> get_s0() const {
      DXTBX_ASSERT(wavelength_ != 0.0);
      return direction_ * 1.0 / wavelength_;
    }

    /** Set the direction and wavelength from s0 */
    void set_s0(vec3<double> s0) {
      DXTBX_ASSERT(s0.length() > 0);
      direction_ = s0.normalize();
      wavelength_ = 1.0 / s0.length();
    }

    /** Check wavlength and direction are (almost) same */
    bool operator==(const Beam &beam) {
      double eps = 1.0e-6;
      double d_direction =  std::abs(angle_safe(direction_, beam.direction_));
      double d_wavelength = std::abs(wavelength_ - beam.wavelength_);
      return (d_direction <= eps && d_wavelength <= eps);
    }

    /** Check wavelength and direction are not (almost) equal. */
    bool operator!=(const Beam &beam) {
      return !(*this == beam);
    }

    friend std::ostream& operator<<(std::ostream &os, const Beam &b);

  private:
    
    double wavelength_;
    vec3 <double> direction_;
  };

  /** Print beam information */
  inline
  std::ostream& operator<<(std::ostream &os, const Beam &b) {
    os << "Beam:\n";
    os << "    wavelength: " << b.get_wavelength() << "\n";
    os << "    direction : " << b.get_direction().const_ref() << "\n";
    return os;
  }

}} // namespace dxtbx::model

#endif // DXTBX_MODEL_BEAM_H
