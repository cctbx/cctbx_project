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

#ifndef DXTBX_MODEL_ATTENUATION_COEFFICIENT_H
#define DXTBX_MODEL_ATTENUATION_COEFFICIENT_H

/** The maximum number of table items */
#define XRAY_MASS_COEFF_TABLE_SIZE 93
#define XRAY_MASS_COEFF_TABLE_MAX_SIZE 66

#include <cmath>
#include <dxtbx/error.h>
#include <scitbx/constants.h>

namespace dxtbx { namespace model {

  using std::log;
  using std::exp;
  using scitbx::constants::factor_ev_angstrom;

  /**
   * The table entries
   */
  struct XrayMassCoeffTableData {
    std::size_t num;
    double energy[XRAY_MASS_COEFF_TABLE_MAX_SIZE];
    double mu_rho[XRAY_MASS_COEFF_TABLE_MAX_SIZE];
    double mu_en_rho[XRAY_MASS_COEFF_TABLE_MAX_SIZE];
  };

  // The array of tables
  extern
  const XrayMassCoeffTableData
    XRAY_MASS_COEFF_TABLE[XRAY_MASS_COEFF_TABLE_SIZE];

  /**
   * A class to access mu_rho and mu_en_rho at different energies
   */
  class XrayMassCoeffTable {
  public:

    /** Construct using pointer to table data */
    XrayMassCoeffTable(std::size_t z) {
      DXTBX_ASSERT(z > 0 && z < XRAY_MASS_COEFF_TABLE_SIZE);
      data_ = &XRAY_MASS_COEFF_TABLE[z];
    }

    /** @return Elements in table */
    std::size_t size() const {
      return data_->num;
    }

    /** @return energy at given index */
    double energy(std::size_t index) const {
      DXTBX_ASSERT(index < size());
      return data_->energy[index];
    }

    /** @return mu_rho at given index */
    double mu_rho(std::size_t index) const {
      DXTBX_ASSERT(index < size());
      return data_->mu_rho[index];
    }

    /** @return mu_en_rho at given index */
    double mu_en_rho(std::size_t index) const {
      DXTBX_ASSERT(index < size());
      return data_->mu_en_rho[index];
    }

    /** @return mu_rho at ev */
    double mu_rho_at_ev(double en) {
      std::size_t index = find_energy_index(en);
      double x0 = log(energy(index));
      double x1 = log(energy(index + 1));
      double y0 = log(mu_rho(index));
      double y1 = log(mu_rho(index + 1));
      return exp(y0 + (y1 - y0) * (log(en) - x0) / (x1 - x0));
    }

    /** @return mu_rho at kev */
    double mu_rho_at_kev(double energy) {
      return mu_rho_at_ev(energy * 1000);
    }

    /** @return mu_rho at wavelength */
    double mu_rho_at_angstrom(double wavelength) {
      return mu_rho_at_ev(factor_ev_angstrom / wavelength);
    }

    /** @return mu_en_rho at ev */
    double mu_en_rho_at_ev(double en) {
      std::size_t index = find_energy_index(en);
      double x0 = log(energy(index));
      double x1 = log(energy(index + 1));
      double y0 = log(mu_en_rho(index));
      double y1 = log(mu_en_rho(index + 1));
      return exp(y0 + (y1 - y0) * (log(en) - x0) / (x1 - x0));
    }

    /** @return mu_en_rho at kev */
    double mu_en_rho_at_kev(double energy) {
      return mu_en_rho_at_ev(energy * 1000);
    }

    /** @return mu_en_rho at wavelength */
    double mu_en_rho_at_angstrom(double wavelength) {
      return mu_en_rho_at_ev(factor_ev_angstrom / wavelength);
    }

  protected:

    /** Find index of energy */
    std::size_t find_energy_index(double energy) {
      std::size_t index = 0;
      for (; index < data_->num; ++index) {
        if (energy < data_->energy[index]) {
          break;
        }
      }
      DXTBX_ASSERT(index > 0 && index < data_->num && data_->num > 0);
      return index - 1;
    }

    const XrayMassCoeffTableData *data_;
  };

}} // namespace dxtbx::model

#endif // DXTBX_MODEL_ATTENUATION_COEFFICIENT_H
