/*
 * robust_glm.h
 *
 *  Copyright (C) 2015 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef SCITBX_GLMTBX_ROBUST_GLM_H
#define SCITBX_GLMTBX_ROBUST_GLM_H

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/matrix/inversion.h>
#include <scitbx/matrix/multiply.h>
#include <scitbx/glmtbx/family.h>

namespace scitbx { namespace glmtbx {

  using scitbx::matrix::inversion_in_place;
  using scitbx::matrix::transpose_multiply;

  namespace detail {

    /**
     * Helper function to compute the sign of a number
     * @param x The number
     * @return The sign -1, 0, 1
     */
    template <typename T>
    int sign(T x) {
      return x == 0 ? 0 : x < 0 ? -1 : 1;
    }

  }

  /**
   * The huber weight function
   * @param r The pearson residual
   * @param c The tuning constant
   * @return The value
   */
  template <typename T>
  T huber(T r, T c) {
    T absr = std::abs(r);
    if (absr < c) {
      return r;
    }
    return detail::sign(r) * c;
  }

  /**
   * A template struct to compute expectation values needed by the algorithm
   */
  template <typename Family>
  struct expectation;

  /**
   * Specialization for poisson distribution
   */
  template <>
  struct expectation<poisson> {

    double epsi1;
    double epsi2;

    /**
     * Compute the expectation values
     * @param mu The mean
     * @param svar The sqrt of the variance
     * @param c The tuning constant
     */
    expectation(double mu, double svar, double c) {

      // Compute some needed probabilities
      int j1 = (int)std::floor(mu - c*svar);
      int j2 = (int)std::floor(mu + c*svar);
      double p1 = poisson::pdf(mu, j1);   // P(Y  = j1)
      double p2 = poisson::pdf(mu, j2);   // P(Y  = j2)
      double p3 = poisson::cdf(mu, j1);   // P(Y <= j1)
      double p4 = poisson::pdf(mu, j2+1); // P(Y  = j2 + 1)
      double p5 = poisson::cdf(mu, j2+1); // P(Y <= j2 + 1)
      double p6 = 1.0 - p5 + p4;          // P(Y >= j2 + 1)
      double p7 = poisson::pdf(mu, j1-1); // P(Y  = j1 - 1)
      double p8 = poisson::pdf(mu, j2-1); // P(Y  = j2 - 1)
      double p9 = poisson::cdf(mu, j2-1); // P(Y <= j2 - 1)
      double p10 = p9 - p3 + p1;          // P(j1 <= Y <= j2)

      // Compute expected values
      epsi1 = c*(p6-p3) + (mu/svar)*(p1-p2);
      epsi2 = c*(p1+p2) + (mu*mu/(svar*svar*svar))*(p10/mu+p7-p1-p8+p2);
    }
  };


  /**
   * An algorithm to do robust generalized linear model as described in
   * Cantoni and Rochetti (2001) "Robust Inference for Generalized Linear
   * Models"
   */
  template <typename Family>
  class robust_glm {

  public:

    typedef Family family;

    /**
     * Compute the generalized linear model using iteratively reweighted least
     * squares. The input expects a design matrix of size (nobs, ncoef), a list
     * of observations of size (nobs) and a list of initial estimates of size
     * (ncoef).
     * @param X The design matrix
     * @param Y The observations
     * @param B The initial estimate
     * @param c The huber tuning constant
     * @param tolerance The stopping critera
     * @param max_iter The maximum number of iterations
     */
    robust_glm(
        const af::const_ref< double, af::c_grid<2> > &X,
        const af::const_ref< double > &Y,
        const af::const_ref< double > &B,
        double c,
        double tolerance,
        std::size_t max_iter)
        : beta_(B.begin(), B.end()),
          niter_(0),
          error_(0),
          c_(c),
          tolerance_(tolerance),
          max_iter_(max_iter) {
      SCITBX_ASSERT(X.accessor()[0] == Y.size());
      SCITBX_ASSERT(X.accessor()[1] == B.size());
      SCITBX_ASSERT(X.accessor()[0] > 0);
      SCITBX_ASSERT(X.accessor()[1] > 0);
      SCITBX_ASSERT(c > 0);
      SCITBX_ASSERT(tolerance > 0);
      SCITBX_ASSERT(max_iter > 0);
      compute(X, Y);
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

    void compute(
        const af::const_ref< double, af::c_grid<2> > &X,
        const af::const_ref< double > &Y) {

      // Number of observations and coefficients
      std::size_t n_obs = X.accessor()[0];
      std::size_t n_cof = X.accessor()[1];

      // Initialize the required matrices and vectors
      af::shared<double> U(n_cof, 0);
      af::versa< double, af::c_grid<2> > BX(af::c_grid<2>(n_obs, n_cof), 0);
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

          // Compute some required values
          double w = 1.0; // TODO Implement different x weights
          double mu  = family::linkinv(eta);
          double var = family::variance(mu);
          double dmu = family::dmu_deta(eta);
          double phi = family::dispersion();
          SCITBX_ASSERT(phi > 0);
          SCITBX_ASSERT(var > 0);
          double svar = std::sqrt(phi * var);
          double res = (Y[i] - mu) / svar;

          // Compute expectation values
          expectation<Family> epsi(mu, svar, c_);

          // Compute the difference between psi and its expected value
          double psi = huber(res, c_);
          double psi_m_epsi = psi - epsi.epsi1;

          // Compute the value of Psi and B_diag for this observation
          double q = psi_m_epsi * w * dmu / svar;
          double b = epsi.epsi2 * w * dmu * dmu / svar;

          // Update the BX = B * X and U matrices
          for (std::size_t j = 0; j < n_cof; ++j) {
            U[j]   += q * X(i,j);
            BX(i,j) = b * X(i,j);
          }
        }

        // Compute the matrix H = X^T B X
        transpose_multiply(
            X.begin(),
            BX.begin(),
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

    af::shared<double> beta_;
    std::size_t niter_;
    double error_;
    double c_;
    double tolerance_;
    std::size_t max_iter_;
  };


}} // namespace scitbx::glmtbx

#endif // SCITBX_GLMTBX_ROBUST_GLM_H
