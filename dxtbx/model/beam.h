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
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/simple_io.h>
#include <scitbx/array_family/simple_tiny_io.h>
#include <dxtbx/error.h>
#include "model_helpers.h"

namespace dxtbx { namespace model {

  using scitbx::vec3;

  /** Base class for beam objects */
  class BeamBase {
  public:
    virtual ~BeamBase() {}

    // Get the direction
    virtual vec3 <double> get_direction() const = 0;
    // Get the wavelength
    virtual double get_wavelength() const = 0;
    // Get the beam divergence
    virtual double get_divergence() const = 0;
    // Get the standard deviation of the beam divergence
    virtual double get_sigma_divergence() const = 0;
    // Set the direction.
    virtual void set_direction(vec3 <double> direction) = 0;
    // Set the wavelength
    virtual void set_wavelength(double wavelength) = 0;
    // Get the wave vector in units of inverse angstroms
    virtual vec3 <double> get_s0() const = 0;
    // Set the direction and wavelength from s0
    virtual void set_s0(vec3<double> s0) = 0;
    // Get the wave vector from source to sample with unit length
    virtual vec3 <double> get_unit_s0() const = 0;
    // Set the direction using the unit_s0 vector
    virtual void set_unit_s0(vec3<double> unit_s0) = 0;
    // Set the beam divergence
    virtual void set_divergence(double divergence) = 0;
    // Set the standard deviation of the beam divergence
    virtual void set_sigma_divergence(double sigma_divergence) = 0;
    // Get the polarization
    virtual vec3 <double> get_polarization_normal() const = 0;
    // Get the polarization fraction
    virtual double get_polarization_fraction() const = 0;
    // Set the polarization plane
    virtual void set_polarization_normal(vec3 <double> polarization_normal) = 0;
    // Set the polarization fraction
    virtual void set_polarization_fraction(double polarization_fraction) = 0;
    // Set the flux
    virtual void set_flux(double flux) = 0;
    // Set the transmission
    virtual void set_transmission(double transmission) = 0;
    // Get the flux
    virtual double get_flux() const = 0;
    // Get the transmission
    virtual double get_transmission() const = 0;
    // @returns the number of scan points
    virtual std::size_t get_num_scan_points() const = 0;
    // Set the s0 vector at scan points
    virtual void set_s0_at_scan_points(const scitbx::af::const_ref< vec3<double> > &s0) = 0;
    // Get the s0 vector at scan points
    virtual scitbx::af::shared< vec3<double> > get_s0_at_scan_points() const = 0;
    // Get the s0 vector at the scan point
    virtual vec3<double> get_s0_at_scan_point(std::size_t index) const = 0;
    // Reset the scan points
    virtual void reset_scan_points() = 0;
    // Check wavelength and direction are (almost) same
    virtual bool operator==(const BeamBase &rhs) const = 0;
    // Check if two models are similar
    virtual bool is_similar_to(
        const BeamBase &rhs,
        double wavelength_tolerance,
        double direction_tolerance,
        double polarization_normal_tolerance,
        double polarization_fraction_tolerance) const = 0;
    // Check wavelength and direction are not (almost) equal.
    virtual bool operator!=(const BeamBase &rhs) const  = 0;
    //  Rotate the beam about an axis
    virtual void rotate_around_origin(vec3<double> axis, double angle) = 0;
  };

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
        polarization_fraction_(0.999),
        flux_(0),
        transmission_(1.0) {}

    /**
     * Initialise all the beam parameters.
     * @param direction The beam direction vector.
     */
    Beam(vec3 <double> s0)
      : divergence_(0.0),
        sigma_divergence_(0.0),
        polarization_normal_(0.0, 1.0, 0.0),
        polarization_fraction_(0.999),
        flux_(0),
        transmission_(1.0) {
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
        polarization_normal_(0.0, 1.0, 0.0),
        polarization_fraction_(0.999),
        flux_(0),
        transmission_(1.0) {
      DXTBX_ASSERT(direction.length() > 0);
      direction_ = direction.normalize();
    }

    /**
     * Initialise all the beam parameters.
     * @param direction The beam direction vector.
     */
    Beam(vec3 <double> s0, double divergence, double sigma_divergence)
      : divergence_(divergence),
        sigma_divergence_(sigma_divergence),
        polarization_normal_(0.0, 1.0, 0.0),
        polarization_fraction_(0.999),
        flux_(0),
        transmission_(1.0) {
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
    Beam(vec3 <double> direction, double wavelength,
         double divergence, double sigma_divergence)
      : wavelength_(wavelength),
        divergence_(divergence),
        sigma_divergence_(sigma_divergence),
        polarization_normal_(0.0, 1.0, 0.0),
        polarization_fraction_(0.999),
        flux_(0),
        transmission_(1.0) {
      DXTBX_ASSERT(direction.length() > 0);
      direction_ = direction.normalize();
    }

    Beam(vec3 <double> direction, double wavelength,
         double divergence, double sigma_divergence,
         vec3<double> polarization_normal,
         double polarization_fraction,
         double flux,
         double transmission)
      : wavelength_(wavelength),
        divergence_(divergence),
        sigma_divergence_(sigma_divergence),
        polarization_normal_(polarization_normal),
        polarization_fraction_(polarization_fraction),
        flux_(0),
        transmission_(1.0) {
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

    /**
     * Set the flux
     */
    void set_flux(double flux) {
      flux_ = flux;
    }

    /**
     * Set the transmission
     */
    void set_transmission(double transmission) {
      transmission_ = transmission;
    }

    /**
     * Get the flux
     */
    double get_flux() const {
      return flux_;
    }

    /**
     * Get the transmission
     */
    double get_transmission() const {
      return transmission_;
    }

    /**
     * @returns the number of scan points
     */
    std::size_t get_num_scan_points() const {
      return s0_at_scan_points_.size();
    }

    /**
     * Set the s0 vector at scan points
     */
    void set_s0_at_scan_points(const scitbx::af::const_ref< vec3<double> > &s0) {
      s0_at_scan_points_ = scitbx::af::shared< vec3<double> >(s0.begin(), s0.end());
    }

    /**
     * Get the s0 vector at scan points
     */
    scitbx::af::shared< vec3<double> > get_s0_at_scan_points() const {
      return s0_at_scan_points_;
    }

    /**
     * Get the s0 vector at the scan point
     */
    vec3<double> get_s0_at_scan_point(std::size_t index) const {
      DXTBX_ASSERT(index < s0_at_scan_points_.size());
      return s0_at_scan_points_[index];
    }

    /**
     * Reset the scan points
     */
    void reset_scan_points() {
      s0_at_scan_points_.clear();
    }

    /** Check wavlength and direction are (almost) same */
    bool operator==(const BeamBase &rhs) const {
      double eps = 1.0e-6;
      return std::abs(angle_safe(direction_, rhs.get_direction())) <= eps
          && std::abs(wavelength_ - rhs.get_wavelength()) <= eps
          && std::abs(divergence_ - rhs.get_divergence()) <= eps
          && std::abs(sigma_divergence_ - rhs.get_sigma_divergence()) <= eps
          && std::abs(angle_safe(polarization_normal_, rhs.get_polarization_normal())) <= eps
          && std::abs(polarization_fraction_ - rhs.get_polarization_fraction()) <= eps;
    }

    /**
     * Check if two models are similar
     */
    bool is_similar_to(
        const BeamBase &rhs,
        double wavelength_tolerance,
        double direction_tolerance,
        double polarization_normal_tolerance,
        double polarization_fraction_tolerance) const {
      return std::abs(angle_safe(direction_, rhs.get_direction())) <= direction_tolerance
          && std::abs(wavelength_ - rhs.get_wavelength()) <= wavelength_tolerance
          && std::abs(angle_safe(polarization_normal_, rhs.get_polarization_normal())) <= polarization_normal_tolerance
          && std::abs(polarization_fraction_ - rhs.get_polarization_fraction()) <= polarization_fraction_tolerance;
    }

    /** Check wavelength and direction are not (almost) equal. */
    bool operator!=(const BeamBase &rhs) const {
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
    double flux_;
    double transmission_;
    scitbx::af::shared< vec3<double> > s0_at_scan_points_;
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
