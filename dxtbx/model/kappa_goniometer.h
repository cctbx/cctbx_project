/*
 * kappa_goniometer.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DXTBX_MODEL_KAPPA_GONIOMETER_H
#define DXTBX_MODEL_KAPPA_GONIOMETER_H

#include <iostream>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/math/r3_rotation.h>
#include <scitbx/constants.h>
#include <dxtbx/error.h>
#include "goniometer.h"

namespace dxtbx { namespace model {

  using scitbx::vec3;
  using scitbx::mat3;
  using scitbx::constants::pi;
  using scitbx::math::r3_rotation::axis_and_angle_as_matrix;

  /**
   * A class representing a kappa goniometer where omega is the primary axis
   * (i,e. aligned with X in the CBF coordinate frame) and has the kappa arm
   * with angle alpha attached to it, aligned with -z, +y, +z or -y at
   * omega = 0, that being the direction, which in turn has phi fixed to it
   * which should initially be coincident with omega. We also need to know
   * which axis is being used for the scan i.e. phi or omega. All angles
   * should be given in degrees. This will work by first constructing the
   * rotation axes and then composing them to the scan axis and fixed
   * component of the rotation.
   */
  class KappaGoniometer : public Goniometer {
  public:

    /** Enumeration of directions */
    enum Direction {
      NoDirection,
      PlusY,
      PlusZ,
      MinusY,
      MinusZ
    };

    /** Enumeration of scan axes */
    enum ScanAxis {
      NoAxis,
      Phi,
      Omega
    };

  public:

    /** Default constructor */
    KappaGoniometer()
      : alpha_(0.0),
        omega_(0.0),
        kappa_(0.0),
        phi_(0.0),
        direction_(NoDirection),
        scan_axis_(NoAxis),
        omega_axis_(0.0, 0.0, 0.0),
        phi_axis_(0.0, 0.0, 0.0),
        kappa_axis_(0.0, 0.0, 0.0) {}

    /**
     * Initialise the goniometer.
     * @param alpha The alpha angle
     * @param omega The omega angle
     * @param kappa The kappa angle
     * @param phi The phi angle
     * @param direction The direction
     * @param scan_axis The scan axis
     */
    KappaGoniometer(double alpha, double omega, double kappa, double phi,
                    Direction direction, ScanAxis scan_axis)
      : alpha_(alpha),
        omega_(omega),
        kappa_(kappa),
        phi_(phi),
        direction_(direction),
        scan_axis_(scan_axis),
        omega_axis_(1.0, 0.0, 0.0),
        phi_axis_(1.0, 0.0, 0.0),
        kappa_axis_(calculate_kappa())
    {
      rotation_axis_ = calculate_rotation_axis();
      fixed_rotation_ = calculate_fixed_rotation();
    }

    /** Virtual destructor */
    virtual ~KappaGoniometer() {}

    /** Get the alpha angle */
    double get_alpha_angle() const {
      return alpha_;
    }

    /** Get the omega angle */
    double get_omega_angle() const {
      return omega_;
    }

    /** Get the kappa angle */
    double get_kappa_angle() const {
      return kappa_;
    }

    /** Get the phi angle */
    double get_phi_angle() const {
      return phi_;
    }

    /** Get the direction */
    Direction get_direction() const {
      return direction_;
    }

    /** Get the scan axis */
    ScanAxis get_scan_axis() const {
      return scan_axis_;
    }

    /** Get the omega axis */
    vec3 <double> get_omega_axis() const {
      return omega_axis_;
    }

    /** Get the phi axis */
    vec3 <double> get_phi_axis() const {
      return phi_axis_;
    }

    /** Get the kappa axis */
    vec3 <double> get_kappa_axis() const {
      return kappa_axis_;
    }

    friend std::ostream& operator<<(std::ostream& os, const KappaGoniometer &g);

  protected:

    /** Calculate the kappa axis */
    vec3 <double> calculate_kappa() {
      double c = cos(alpha_ * pi / 180.0);
      double s = sin(alpha_ * pi / 180.0);
      if (direction_ == PlusY) {
        return vec3 <double> (c, s, 0.0);
      } else if (direction_ == PlusZ) {
        return vec3 <double> (c, 0.0, s);
      } else if (direction_ == MinusY) {
        return vec3 <double> (c, -s, 0.0);
      } else if (direction_ == MinusZ) {
        return vec3 <double> (c, 0.0, -s);
      } else {
        throw DXTBX_ERROR("Invalid direction");
      }
      return vec3 <double> (0.0, 0.0, 0.0);
    }

    /** Calculate the rotation axis */
    vec3 <double> calculate_rotation_axis() {
      if (scan_axis_ == Omega) {
        return omega_axis_;
      } else if (scan_axis_ == Phi) {
        mat3 <double> O = axis_and_angle_as_matrix(omega_axis_, omega_, true);
        mat3 <double> K = axis_and_angle_as_matrix(kappa_axis_, kappa_, true);
        return O * K * phi_axis_;
      } else {
        throw DXTBX_ERROR("Invalid scan axis");
      }
      return vec3 <double> (0.0, 0.0, 0.0);
    }

    /** Calculate the fixed rotation */
    mat3 <double> calculate_fixed_rotation() {
      if (scan_axis_ == Omega) {
        mat3 <double> K = axis_and_angle_as_matrix(kappa_axis_, kappa_, true);
        mat3 <double> P = axis_and_angle_as_matrix(phi_axis_, phi_, true);
        return K * P;
      } else if (scan_axis_ == Phi) {
        return mat3 <double> (
          1.0, 0.0, 0.0,
          0.0, 1.0, 0.0,
          0.0, 0.0, 1.0);
      } else {
        throw DXTBX_ERROR("Invalid scan axis");
      }
      return mat3 <double> (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    }

    double alpha_;
    double omega_;
    double kappa_;
    double phi_;
    Direction direction_;
    ScanAxis scan_axis_;
    vec3 <double> omega_axis_;
    vec3 <double> phi_axis_;
    vec3 <double> kappa_axis_;
  };

  /** Print direction enum */
  inline
  std::ostream& operator<<(std::ostream &os, KappaGoniometer::Direction d) {
    if (d == KappaGoniometer::PlusY) {
      os << "+y";
    } else if (d == KappaGoniometer::PlusZ) {
      os << "+z";
    } else if (d == KappaGoniometer::MinusY) {
      os << "-y";
    } else if (d == KappaGoniometer::MinusZ) {
      os << "-z";
    } else {
      os << "none";
    }
    return os;
  }

  /** Print scan axis enum */
  inline
  std::ostream& operator<<(std::ostream &os, KappaGoniometer::ScanAxis a) {
    if (a == KappaGoniometer::Omega) {
      os << "omega";
    } else if (a == KappaGoniometer::Phi) {
      os << "phi";
    } else {
      os << "none";
    }
    return os;
  }

  /** Print kappa goniometer info */
  inline
  std::ostream& operator<<(std::ostream& os, const KappaGoniometer &g) {
    os << "Goniometer:\n";
    os << "    Rotation axis:  " << g.get_rotation_axis().const_ref() << "\n";
    os << "    Fixed rotation: " << g.get_fixed_rotation().const_ref() << "\n";
    os << "    alpha angle:    " << g.get_alpha_angle() << "\n";
    os << "    omega angle:    " << g.get_omega_angle() << "\n";
    os << "    kappa angle:    " << g.get_kappa_angle() << "\n";
    os << "    phi angle:      " << g.get_phi_angle() << "\n";
    os << "    direction:      " << g.get_direction() << "\n";
    os << "    scan axis:      " << g.get_scan_axis() << "\n";
    os << "    omega axis:     " << g.get_omega_axis().const_ref() << "\n";
    os << "    kappa axis:     " << g.get_kappa_axis().const_ref() << "\n";
    os << "    phi axis:       " << g.get_phi_axis().const_ref() << "\n";
    return os;
  }


}} // namespace dxtbx::model

#endif // DXTBX_MODEL_KAPPA_GONIOMETER_H
