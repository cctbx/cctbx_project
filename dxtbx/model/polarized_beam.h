/*
 * polarized_beam.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DXTBX_MODEL_POLARIZED_BEAM_H
#define DXTBX_MODEL_POLARIZED_BEAM_H

#include <cmath>
#include <iostream>
#include <scitbx/vec3.h>
#include "beam.h"

namespace dxtbx { namespace model {

  using scitbx::vec3;

  /** A class to represent a polarized beam. */
  class PolarizedBeam : public Beam {
  public:
    /** Default constructor: initialise all to zero */
    PolarizedBeam()
      : Beam(),
        polarization_(0.0, 0.0, 0.0),
        polarization_fraction_(0.0) {}

    /**
     * Initialise all the beam parameters.
     * @param wavelength The wavelength of the beam
     * @param direction The beam direction vector.
     * @param polarization The polarization plane of the beam
     * @param polarization_fraction The polarization fraction.
     */
    PolarizedBeam(vec3 <double> direction,
                  double wavelength,
                  vec3 <double> polarization,
                  double polarization_fraction)
      : Beam(direction, wavelength),
        polarization_(polarization),
        polarization_fraction_(polarization_fraction) {}

    /**
     * Initialise all the beam parameters.
     * @param direction The beam direction vector.
     * @param polarization The polarization plane of the beam
     * @param polarization_fraction The polarization fraction.
     */
    PolarizedBeam(vec3 <double> s0,
                  vec3 <double> polarization,
                  double polarization_fraction)
      : Beam(s0),
        polarization_(polarization),
        polarization_fraction_(polarization_fraction) {}

    /** Virtual destructor */
    virtual ~PolarizedBeam() {}

    /** Get the polarization */
    vec3 <double> get_polarization() const {
      return polarization_;
    }

    /** Get the polarization fraction */
    double get_polarization_fraction() const {
      return polarization_fraction_;
    }

    /** Set the polarization plane */
    void set_polarization(vec3 <double> polarization) {
      polarization_ = polarization;
    }

    /** Set the polarization fraction */
    void set_polarization_fraction(double polarization_fraction) {
      polarization_fraction_ = polarization_fraction;
    }

    friend std::ostream& operator<<(std::ostream &os, const PolarizedBeam &b);

  private:
    vec3 <double> polarization_;
    double polarization_fraction_;
  };

  /** Print the beam info */
  inline
  std::ostream& operator<<(std::ostream &os, const PolarizedBeam &b) {
    os << "Beam:\n";
    os << "    wavelength:   " << b.get_wavelength() << "\n";
    os << "    direction :   " << b.get_direction().const_ref() << "\n";
    os << "    Polarization: " << b.get_polarization().const_ref() << "\n";
    os << "    Pn fraction:  " << b.get_polarization_fraction() << "\n";
    return os;
  }

}} // namespace dxtbx::model

#endif // DXTBX_MODEL_POLARIZED_BEAM_H
