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
        direction_(0.0, 0.0, 0.0),
        divergence_(0.0),
        sigma_divergence_(0.0),
        polarization_normal_(0.0, 1.0, 0.0),
        polarization_fraction_(0.999) {}

    /**
     * Initialise all the beam parameters.
     * @param direction The beam direction vector.
     */
    Beam(vec3 <double> s0)
      : divergence_(0.0),
        sigma_divergence_(0.0),
        polarization_normal_(0.0, 1.0, 0.0),
        polarization_fraction_(0.999) {
      DXTBX_ASSERT(s0.length() > 0);
      wavelength_ = 1.0 / s0.length();
      direction_ = -s0.normalize();
    }

    /**
     * Initialise all the beam parameters. Normalize the direction vector
     * and give it the length of 1.0 / wavelength
     * @param wavelength The wavelength of the beam
     * @param direction The beam direction vector.
     */
    Beam(vec3 <double> direction, double wavelength)
      : wavelength_(wavelength),
        divergence_(0.0),
        sigma_divergence_(0.0),
        polarization_fraction_(0.999) {
      DXTBX_ASSERT(direction.length() > 0);
      direction_ = direction.normalize();
      polarization_normal_ = direction_ * vec3<double>(1, 0, 0);
    }

    /**
     * Initialise all the beam parameters.
     * @param direction The beam direction vector.
     */
    Beam(vec3 <double> s0, double divergence, double sigma_divergence)
      : divergence_(divergence),
        sigma_divergence_(sigma_divergence),
        polarization_fraction_(0.999) {
      DXTBX_ASSERT(s0.length() > 0);
      wavelength_ = 1.0 / s0.length();
      direction_ = -s0.normalize();
      polarization_normal_ = direction_ * vec3<double>(1, 0, 0);
    }

    /**
     * Initialise all the beam parameters. Normalize the direction vector
     * and give it the length of 1.0 / wavelength
     * @param wavelength The wavelength of the beam
     * @param direction The beam direction vector.
     */
    Beam(vec3 <double> direction, double wavelength,
         double divergence, double sigma_divergence)
      : wavelength_(wavelength),
        divergence_(divergence),
        sigma_divergence_(sigma_divergence),
        polarization_fraction_(0.999)  {
      DXTBX_ASSERT(direction.length() > 0);
      direction_ = direction.normalize();
      polarization_normal_ = direction_ * vec3<double>(1, 0, 0);
    }

    Beam(vec3 <double> direction, double wavelength,
         double divergence, double sigma_divergence,
         vec3<double> polarization_normal,
         double polarization_fraction)
      : wavelength_(wavelength),
        divergence_(divergence),
        sigma_divergence_(sigma_divergence),
        polarization_normal_(polarization_normal),
        polarization_fraction_(polarization_fraction)  {
      direction_ = direction.normalize();
      const double EPS = 1e-7;
      DXTBX_ASSERT(std::abs(direction_ * polarization_normal_) < EPS);
      DXTBX_ASSERT(direction.length() > 0);
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

    /** Get the beam divergence */
    double get_divergence() const {
      return divergence_;
    }

    /** Get the standard deviation of the beam divergence */
    double get_sigma_divergence() const {
      return sigma_divergence_;
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
      return -direction_ * 1.0 / wavelength_;
    }

    /** Set the direction and wavelength from s0 */
    void set_s0(vec3<double> s0) {
      DXTBX_ASSERT(s0.length() > 0);
      direction_ = -s0.normalize();
      wavelength_ = 1.0 / s0.length();
    }

    /** Get the wave vector from source to sample with unit length */
    vec3 <double> get_unit_s0() const {
      return -direction_;
    }

    /** Set the direction using the unit_s0 vector */
    void set_unit_s0(vec3<double> unit_s0) {
      DXTBX_ASSERT(unit_s0.length() > 0);
      direction_ = -(unit_s0.normalize());
    }

    /** Set the beam divergence */
    void set_divergence(double divergence) {
      divergence_ = divergence;
    }

    /** Set the standard deviation of the beam divergence */
    void set_sigma_divergence(double sigma_divergence) {
      sigma_divergence_ = sigma_divergence;
    }

    /** Get the polarization */
    vec3 <double> get_polarization_normal() const {
      return polarization_normal_;
    }

    /** Get the polarization fraction */
    double get_polarization_fraction() const {
      return polarization_fraction_;
    }

    /** Set the polarization plane */
    void set_polarization_normal(vec3 <double> polarization_normal) {
      polarization_normal_ = polarization_normal;
    }

    /** Set the polarization fraction */
    void set_polarization_fraction(double polarization_fraction) {
      polarization_fraction_ = polarization_fraction;
    }

    /** Check wavlength and direction are (almost) same */
    bool operator==(const Beam &rhs) {
      double eps = 1.0e-6;
      return std::abs(angle_safe(direction_, rhs.direction_)) <= eps
          && std::abs(wavelength_ - rhs.wavelength_) <= eps
          && std::abs(divergence_ - rhs.divergence_) <= eps
          && std::abs(sigma_divergence_ - rhs.sigma_divergence_) <= eps
          && std::abs(angle_safe(polarization_normal_, rhs.polarization_normal_)) <= eps
          && std::abs(polarization_fraction_ - rhs.polarization_fraction_) <= eps;
    }

    /** Check wavelength and direction are not (almost) equal. */
    bool operator!=(const Beam &rhs) {
      return !(*this == rhs);
    }

    /** Rotate the beam about an axis */
    void rotate_around_origin(vec3<double> axis, double angle) {
      const double EPS = 1e-7;
      DXTBX_ASSERT(std::abs(direction_ * polarization_normal_) < EPS);
      direction_ = direction_.rotate_around_origin(axis, angle);
      polarization_normal_ = polarization_normal_.rotate_around_origin(axis, angle);
    }

    friend std::ostream& operator<<(std::ostream &os, const Beam &b);

  private:

    double wavelength_;
    vec3 <double> direction_;
    double divergence_;
    double sigma_divergence_;
    vec3 <double> polarization_normal_;
    double polarization_fraction_;
  };

  /** Print beam information */
  inline
  std::ostream& operator<<(std::ostream &os, const Beam &b) {
    os << "Beam:\n";
    os << "    wavelength: " << b.get_wavelength() << "\n";
    os << "    sample to source direction : " << b.get_direction().const_ref() << "\n";
    os << "    divergence: " << b.get_divergence() << "\n";
    os << "    sigma divergence: " << b.get_sigma_divergence() << "\n";
    os << "    polarization normal: " <<
        b.get_polarization_normal().const_ref() << "\n";
    os << "    polarization fraction: " <<
        b.get_polarization_fraction() << "\n";
    return os;
  }

}} // namespace dxtbx::model

#endif // DXTBX_MODEL_BEAM_H
