#ifndef SCITBX_MINPACK_LEVENBERG_MARQUARDT_H
#define SCITBX_MINPACK_LEVENBERG_MARQUARDT_H

#include <scitbx/minpack/raw.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/error.h>
#include <vector>

namespace scitbx { namespace minpack {

  //! C++ wrapper for MINPACK lmder.
  class levenberg_marquardt
  {
    public:
      //! Default constructor. Some data members are not initialized!
      levenberg_marquardt() {}

      //! Initialization of the minimization.
      levenberg_marquardt(
        int m_,
        af::shared<double> x_,
        double ftol_=-1,
        double xtol_=-1,
        double gtol_=0,
        int maxfev_=0,
        double factor_=1.0e2,
        bool call_back_after_iteration=false)
      :
        m(m_),
        x(x_),
        ftol(ftol_),
        xtol(xtol_),
        gtol(gtol_),
        maxfev(maxfev_),
        factor(factor_),
        fvec(m_, 0),
        fjac(m_*x_.size(), 0),
        ipvt(x_.size(), 0),
        workspace(5*x_.size()+m_, 0),
        info(0),
        nfev(0),
        njev(0),
        minimizer(call_back_after_iteration)
      {
        static const double default_tol = std::sqrt(raw::dpmpar(1));
        if (ftol < 0) ftol = default_tol;
        if (xtol < 0) xtol = ftol;
        run();
        SCITBX_ASSERT(minimizer.iflag == 1);
      }

      bool
      has_terminated() const
      {
        return minimizer.iflag == 0;
      }

      bool
      requests_fvec() const
      {
        return minimizer.iflag == 1 || minimizer.iflag == 3;
      }

      bool
      requests_fjac() const
      {
        return minimizer.iflag == 2;
      }

      bool
      calls_back_after_iteration() const
      {
        return minimizer.iflag == 4;
      }

      void
      process_fvec(af::const_ref<double> const& fvec)
      {
        SCITBX_ASSERT(requests_fvec());
        SCITBX_ASSERT(fvec.size() == m);
        if (minimizer.iflag == 3) {
          std::copy(x_buffer.begin(), x_buffer.end(), x.begin());
          std::copy(fvec.begin(), fvec.end(),
            (workspace.size() ? &*workspace.begin()+5*x.size() : 0));
        }
        else {
          std::copy(fvec.begin(), fvec.end(), this->fvec.begin());
        }
        run();
      }

      void
      process_fjac(af::const_ref<double> const& fjac)
      {
        SCITBX_ASSERT(requests_fjac());
        SCITBX_ASSERT(fjac.size() == m * x.size());
        std::copy(fjac.begin(), fjac.end(), this->fjac.begin());
        run();
        SCITBX_ASSERT(has_terminated()
                   || requests_fvec()
                   || calls_back_after_iteration());
      }

      void
      continue_after_call_back_after_iteration()
      {
        SCITBX_ASSERT(calls_back_after_iteration());
        run();
      }

      //! Not available from Python.
      void
      run()
      {
        int mode = 1;
        int x_size = static_cast<int>(x.size());
        double* ws0 = (workspace.size() ? &*workspace.begin() : 0);
        minimizer.run(
          m,
          x_size,
          raw::ref1<double>(x.begin(), x_size),
          raw::ref1<double>(fvec.begin(), fvec.size()),
          raw::ref2<double>(fjac.begin(), m, x_size),
          m,
          ftol,
          xtol,
          gtol,
          maxfev,
          raw::ref1<double>(ws0, x_size),
          mode,
          factor,
          info,
          nfev,
          njev,
          raw::ref1<int>((ipvt.size() ? &*ipvt.begin() : 0), ipvt.size()),
          raw::ref1<double>(ws0+1*x_size, x_size),
          raw::ref1<double>(ws0+2*x_size, x_size),
          raw::ref1<double>(ws0+3*x_size, x_size),
          raw::ref1<double>(ws0+4*x_size, x_size),
          raw::ref1<double>(ws0+5*x_size, m));
        if (minimizer.iflag == 3) {
          x_buffer.assign(x.begin(), x.end());
          std::copy(ws0+3*x_size, ws0+4*x_size, x.begin());
        }
      }

      int m;
      af::shared<double> x;
      double ftol;
      double xtol;
      double gtol;
      int maxfev;
      double factor;
      af::shared<double> fvec;
      af::shared<double> fjac;
      std::vector<int> ipvt;
      std::vector<double> workspace;
      std::vector<double> x_buffer;
      int info;
      int nfev;
      int njev;
      raw::lmder minimizer;
  };

}} // namespace scitbx::minpack

#endif // SCITBX_MINPACK_LEVENBERG_MARQUARDT_H
