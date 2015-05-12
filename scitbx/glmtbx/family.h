/*
 * family.h
 *
 *  Copyright (C) 2015 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef SCITBX_GLMTBX_FAMILY_H
#define SCITBX_GLMTBX_FAMILY_H

#include <boost/math/distributions/poisson.hpp>
#include <cmath>

namespace scitbx { namespace glmtbx {

  /**
   * Define some functions needed to do a poisson GLM
   */
  class poisson {
  public:

    /**
     * The link function (log)
     * @param mu The poisson parameters
     * @return log(mu)
     */
    static
    double link(double mu) {
      return std::log(mu);
    }

    /**
     * The inverse link function
     * @param eta The value of the link
     * @return The poisson parameter
     */
    static
    double linkinv(double eta) {
      return std::exp(eta);
    }

    /**
     * The variance as a function of mu
     * @param mu The poisson parameter
     * @return The variance
     */
    static
    double variance(double mu) {
      return mu;
    }

    /**
     * The derivative of mu wrt eta
     * @param eta
     * @return The value of the derivative
     */
    static
    double dmu_deta(double eta) {
      return std::exp(eta);
    }

    /**
     * The derivative of eta wrt mu
     * @param mu
     * @return The value of the derivative
     */
    static
    double deta_dmu(double mu) {
      SCITBX_ASSERT(mu > 0);
      return 1.0 / mu;
    }

    /**
     * The dispersion
     * @return The dispersion
     */
    static
    double dispersion() {
      return 1.0;
    }

    /**
     * The pdf at the given value
     * @param mean The poisson mean
     * @param value The value at which compute the pdf
     * @return The pdf value
     */
    static
    double pdf(double mean, double value) {
      if (mean == 0) {
        return 0.0;
      }
      if (value == 0) {
        return std::exp(-mean);
      }
      if (value < 0) {
        return 0;
      }
      return std::exp(value * std::log(mean) - mean - boost::math::lgamma(value+1));
    }

    /**
     * The cdf at the given value
     * @param mean The poisson mean
     * @param value The value at which compute the cdf
     * @return The cdf value
     */
    static
    double cdf(double mean, double value) {
      if (mean == 0) {
        return 0.0;
      }
      if (value == 0) {
        return std::exp(-mean);
      }
      if (value < 0) {
        return 0;
      }
      return boost::math::gamma_q(std::floor(value+1), mean);
    }

  };



}} // namespace scitbx::glmtbx


#endif // SCITBX_GLMTBX_FAMILY_H
