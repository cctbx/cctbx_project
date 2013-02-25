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

#include <cmath>
#include <scitbx/vec3.h>

namespace dxtbx { namespace model {

  using scitbx::vec3;

  /** Base class for beam objects */
  class BeamBase {};

  /** A class to represent a simple beam. */
  class Beam : public BeamBase {
  public:
    /** Default constructor: initialise all to zero */
    Beam()
      : direction_(0.0, 0.0, 0.0),
        wavelength_(0.0) {}

    /**
     * Initialise all the beam parameters.
     * @param direction The beam direction vector.
     */
    Beam(vec3 <double> direction)
      : direction_(direction),
        wavelength_(1.0 / direction.length()) {}

    /**
     * Initialise all the beam parameters. Normalize the direction vector
     * and give it the length of 1.0 / wavelength
     * @param wavelength The wavelength of the beam
     * @param direction The beam direction vector.
     */
    Beam(vec3 <double> direction, double wavelength)
      : direction_(direction.normalize() / wavelength),
        wavelength_(wavelength) {}

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
      direction_ = direction;
      wavelength_ = 1.0 / direction.length();
    }

    /** Check wavlength and direction are (almost) same */
    bool operator==(const Beam &beam) {
      double eps = 1.0e-6;
      double d_direction =  std::abs(direction_.angle(beam.direction_));
      double d_wavelength = std::abs(wavelength_ - beam.wavelength_);
      return (d_direction <= eps && d_wavelength <= eps);
    }

    /** Check wavelength and direction are not (almost) equal. */
    bool operator!=(const Beam &beam) {
      return !(*this == beam);
    }

  private:
    vec3 <double> direction_;
    double wavelength_;
  };

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
    PolarizedBeam(vec3 <double> direction,
                  vec3 <double> polarization,
                  double polarization_fraction)
      : Beam(direction),
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

  private:
    vec3 <double> polarization_;
    double polarization_fraction_;
};

}} // namespace dxtbx::model

#endif // DXTBX_MODEL_BEAM_H
