/*
 * glm.h
 *
 *  Copyright (C) 2015 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef SCITBX_GLMTBX_GLM_H
#define SCITBX_GLMTBX_GLM_H

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/matrix/inversion.h>
#include <scitbx/matrix/multiply.h>
#include <scitbx/glmtbx/family.h>

namespace scitbx { namespace glmtbx {

  using scitbx::matrix::inversion_in_place;
  using scitbx::matrix::transpose_multiply;

  /**
   * Compute the generalized linear model using the given family
   */
  template <typename Family>
  class glm {

  public:

    typedef Family family;

    /**
     * Compute the generalized linear model using iteratively reweighted least
     * squares. The input expects a design matrix of size (nobs, ncoef), a list
     * of observations of size (nobs) and a list of initial estimates of size
     * (ncoef) and a list of weights of size (nobs).
     */
    glm(
        const af::const_ref< double, af::c_grid<2> > &X,
        const af::const_ref< double > &Y,
        const af::const_ref< double > &B,
        const af::const_ref< double > &P,
        double tolerance,
        std::size_t max_iter)
        : beta_(B.begin(), B.end()),
          tolerance_(tolerance),
          max_iter_(max_iter),
          error_(0),
          niter_(0) {
      SCITBX_ASSERT(X.accessor()[0] == Y.size());
      SCITBX_ASSERT(X.accessor()[1] == B.size());
      SCITBX_ASSERT(X.accessor()[0] == P.size());
      SCITBX_ASSERT(X.accessor()[0] > 0);
      SCITBX_ASSERT(X.accessor()[1] > 0);
      SCITBX_ASSERT(tolerance > 0);
      SCITBX_ASSERT(max_iter > 0);

      // Number of observations and coefficients
      std::size_t n_obs = X.accessor()[0];
      std::size_t n_cof = X.accessor()[1];

      // Initialize some stuff
      af::shared<double> U(n_cof, 0);
      af::versa< double, af::c_grid<2> > WX(af::c_grid<2>(n_obs, n_cof), 0);
      af::versa< double, af::c_grid<2> > H(af::c_grid<2>(n_cof, n_cof), 0);

      // Loop until we reach the maximum number of iterations
      for (niter_ = 0; niter_ < max_iter_; ++niter_) {

        // Initialize the sum to zero
        for (std::size_t j = 0; j < n_cof; ++j) {
          U[j] = 0.0;
        }

        // Build the matrices from the observations
        for (std::size_t i = 0; i < n_obs; ++i) {

          // Compute the values for eta
          double eta = 0.0;
          for (std::size_t j = 0; j < n_cof; ++j) {
            eta += X(i,j) * beta_[j];
          }

          // Compute some required stuff
          double mu  = family::linkinv(eta);
          double var = family::variance(mu);
          double phi = family::dispersion();
          double deta_dmu = family::deta_dmu(mu);
          SCITBX_ASSERT(deta_dmu > 0);
          SCITBX_ASSERT(phi > 0);
          SCITBX_ASSERT(var > 0);

          // Construct the weight matrix and residual vector
          double w = P[i] / (var*deta_dmu*deta_dmu);
          double z = (Y[i] - mu) * deta_dmu;

          // Update the WX = B * X and U matrices
          for (std::size_t j = 0; j < n_cof; ++j) {
            U[j]   += X(i,j) * w * z;
            WX(i,j) = w * X(i,j);
          }
        }

        // Compute the matrix H = X^T W X
        transpose_multiply(
            X.begin(),
            WX.begin(),
            n_obs,
            n_cof,
            n_cof,
            H.begin());

        // Compute delta = H^-1 U
        inversion_in_place(
            H.begin(),
            n_cof,
            U.begin(),
            1);

        // Compute the relative error in the parameters and update
        double sum_delta_sq = 0.0;
        double sum_beta_sq = 0.0;
        for (std::size_t j = 0; j < n_cof; ++j) {
          sum_delta_sq += U[j]*U[j];
          sum_beta_sq += beta_[j]*beta_[j];
          beta_[j] += U[j];
        }

        // If error is within tolerance then break
        error_ = std::sqrt(sum_delta_sq / std::max(1e-10, sum_beta_sq));
        if (error_ < tolerance_) {
          break;
        }
      }
    }

    /**
     * @returns The parameters
     */
    af::shared<double> parameters() const {
      return beta_;
    }

    /**
     * @returns The number of iterations
     */
    std::size_t niter() const {
      return niter_;
    }

    /**
     * @returns The reletive error at the last iteration
     */
    double error() const {
      return error_;
    }

    /**
     * @returns Did the algorithm converge
     */
    bool converged() const {
      return niter_ < max_iter_;
    }

    /**
     * Compute the values of mu, at X given the computed parameters
     * @param X The design matrix
     * @return The values
     */
    af::shared<double> mu(
        const af::const_ref< double, af::c_grid<2> > &X) const {
      SCITBX_ASSERT(X.accessor()[1] == beta_.size());
      af::shared<double> result(X.accessor()[0]);
      for (std::size_t i = 0; i < result.size(); ++i) {
        double eta = 0.0;
        for (std::size_t j = 0; j < beta_.size(); ++j) {
          eta += X(i,j) * beta_[j];
        }
        result[i] = family::linkinv(eta);
      }
      return result;
    }

  private:

    af::shared<double> beta_;
    std::size_t niter_;
    double error_;
    double tolerance_;
    std::size_t max_iter_;

  };

}} // namespace scitbx::glmtbx


#endif // SCITBX_GLMTBX_GLM_H
