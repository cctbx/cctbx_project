/*
 * attenuation_coefficient.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef CCTBX_ELTBX_ATTENUATION_COEFFICIENT_H
#define CCTBX_ELTBX_ATTENUATION_COEFFICIENT_H

/** The maximum number of table items */
#define XRAY_MASS_COEFF_TABLE_SIZE 95
#define XRAY_MASS_COEFF_TABLE_MAX_SIZE 66

#include <cmath>
#include <cctbx/error.h>
#include <scitbx/constants.h>
#include <scitbx/array_family/flex_types.h>

namespace cctbx { namespace eltbx { namespace attenuation_coefficient {

  using std::log;
  using std::exp;
  using scitbx::constants::factor_ev_angstrom;
  using scitbx::af::flex_double;

  /**
   * The table entries
   */
  struct table_data {
    std::size_t num;
    double energy[XRAY_MASS_COEFF_TABLE_MAX_SIZE];
    double mu_rho[XRAY_MASS_COEFF_TABLE_MAX_SIZE];
    double mu_en_rho[XRAY_MASS_COEFF_TABLE_MAX_SIZE];
  };

  // The array of densities
  extern
  const double DENSITY[XRAY_MASS_COEFF_TABLE_SIZE];

  // The array of tables
  extern
  const table_data XRAY_MASS_COEFF_TABLE[XRAY_MASS_COEFF_TABLE_SIZE];

  /**
   * A class to access mu_rho and mu_en_rho at different energies
   */
  class table {
  public:

    /** Construct using pointer to table data */
    table(std::size_t z) {
      CCTBX_ASSERT(z > 0 && z < XRAY_MASS_COEFF_TABLE_SIZE);
      data_ = &XRAY_MASS_COEFF_TABLE[z];
      density_ = DENSITY[z];
    }

    /** @return Elements in table */
    std::size_t size() const {
      return data_->num;
    }

    /** @return the density in g/cm^3 */
    double density() const {
      return density_;
    }

    /** @return The minimum energy in mev */
    double min_energy() const {
      return energy(0);
    }

    /** @return The maximum energy in mec */
    double max_energy() const {
      return energy(size() - 1);
    }

    /** @return energy at given index in Mev */
    double energy(std::size_t index) const {
      CCTBX_ASSERT(index < size());
      return data_->energy[index];
    }

    /** @return mu_rho at given index in cm^2 / g */
    double mu_rho(std::size_t index) const {
      CCTBX_ASSERT(index < size());
      return data_->mu_rho[index];
    }

    /** @return mu_en_rho at given index in cm^2 / g */
    double mu_en_rho(std::size_t index) const {
      CCTBX_ASSERT(index < size());
      return data_->mu_en_rho[index];
    }

    /** @return The list of energies in Mev */
    flex_double energy() const {
      flex_double result(size());
      for (std::size_t i = 0; i < size(); ++i) {
        result[i] = data_->energy[i];
      }
      return result;
    }

    /** @return The list of mu rho in cm^2 / g */
    flex_double mu_rho() const {
      flex_double result(size());
      for (std::size_t i = 0; i < size(); ++i) {
        result[i] = data_->mu_rho[i];
      }
      return result;
    }

    /** @return The lust of mu_en_rho in cm^2 / g */
    flex_double mu_en_rho() const {
      flex_double result(size());
      for (std::size_t i = 0; i < size(); ++i) {
        result[i] = data_->mu_en_rho[i];
      }
      return result;
    }

    /** @return mu_rho at ev in cm^2 / g */
    double mu_rho_at_ev(double en) const {
      en = en / 1000000;
      std::size_t index = find_energy_index(en);
      double x0 = log(energy(index));
      double x1 = log(energy(index + 1));
      double y0 = log(mu_rho(index));
      double y1 = log(mu_rho(index + 1));
      return exp(y0 + (y1 - y0) * (log(en) - x0) / (x1 - x0));
    }

    /** @return mu_rho at kev in cm^2 / g */
    double mu_rho_at_kev(double energy) const {
      return mu_rho_at_ev(energy * 1000);
    }

    /** @return mu_rho at wavelength in cm^2 / g */
    double mu_rho_at_angstrom(double wavelength) const {
      return mu_rho_at_ev(factor_ev_angstrom / wavelength);
    }

    /** @return mu_en_rho at ev in cm^2 / g */
    double mu_en_rho_at_ev(double en) const {
      en = en / 1000000;
      std::size_t index = find_energy_index(en);
      double x0 = log(energy(index));
      double x1 = log(energy(index + 1));
      double y0 = log(mu_en_rho(index));
      double y1 = log(mu_en_rho(index + 1));
      return exp(y0 + (y1 - y0) * (log(en) - x0) / (x1 - x0));
    }

    /** @return mu_en_rho at kev in cm^2 / g */
    double mu_en_rho_at_kev(double energy) const {
      return mu_en_rho_at_ev(energy * 1000);
    }

    /** @return mu_en_rho at wavelength in cm^2 / g */
    double mu_en_rho_at_angstrom(double wavelength) const {
      return mu_en_rho_at_ev(factor_ev_angstrom / wavelength);
    }

    /** @return mu at ev in cm^-1 */
    double mu_at_ev(double en) const {
      return mu_rho_at_ev(en) * density_;
    }

    /** @return mu at kev in cm^-1 */
    double mu_at_kev(double energy) const {
      return mu_at_ev(energy * 1000);
    }

    /** @return mu at wavelength in cm^-1 */
    double mu_at_angstrom(double wavelength) const {
      return mu_at_ev(factor_ev_angstrom / wavelength);
    }

    /** @return The attenuation length in cm */
    double lambda_at_ev(double energy) const {
      return 1.0 / mu_at_ev(energy);
    }

    /** @return The attenuation length in cm */
    double lambda_at_kev(double energy) const {
      return 1.0 / mu_at_kev(energy);
    }

    /** @return The attenuation length in cm */
    double lambda_at_angstrom(double wavelength) const {
      return 1.0 / mu_at_angstrom(wavelength);
    }

  protected:

    /** Find index of energy */
    std::size_t find_energy_index(double energy) const {
      std::size_t index = 0;
      for (; index < data_->num; ++index) {
        if (energy < data_->energy[index]) {
          break;
        }
      }
      CCTBX_ASSERT(index > 0 && index < data_->num && data_->num > 0);
      return index - 1;
    }

    const table_data *data_;
    double density_;
  };

}}} // namespace cctbx::eltbx::attenuation_coefficient

#endif // CCTBX_ELTBX_ATTENUATION_COEFFICIENT_H
