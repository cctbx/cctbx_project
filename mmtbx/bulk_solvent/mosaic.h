#ifndef MMTBX_MOSAIC_H
#define MMTBX_MOSAIC_H

#include <mmtbx/import_scitbx_af.h>
#include <scitbx/matrix/eigensystem.h>
#include <scitbx/array_family/versa_matrix.h>
#include <boost/python/list.hpp>
#include <scitbx/math/utils.h>
#include <mmtbx/error.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <vector>

namespace mmtbx { namespace mosaic {

typedef std::complex<double> ComplexType;


template <typename FloatType=double>
class alg2_tg
{
public:
  alg2_tg() {}

  alg2_tg(
    boost::python::list      const& F_,
    af::const_ref<FloatType> const& i_obs_)
  :
  dim(boost::python::len(F_)),
  size(i_obs_.size()),
  target_(0),
  gradient_(boost::python::len(F_)),
  F_conj(boost::python::len(F_)),
  F(boost::python::len(F_)),
  i_obs(i_obs_.size()),
  i_obs_sum(0)
  {
   for(std::size_t i=0; i < size; i++) {
     i_obs[i] = i_obs_[i];
     i_obs_sum += i_obs_[i];
   }
   // Extract Fs
   for(std::size_t i=0; i < dim; i++) {
     boost::python::extract<af::const_ref<ComplexType> > elem_proxy_1(F_[i]);
     af::const_ref<ComplexType> fm = elem_proxy_1();
     F[i] = fm;
     MMTBX_ASSERT(fm.size() == i_obs.size());
   }
   // Pre-compute complex conj
   for(std::size_t i=0; i < dim; i++) {
     af::shared<ComplexType> f(size);
     for(std::size_t j=0; j < size; j++) {
       f[j]=std::conj(F[i][j]);
     }
     F_conj[i] = f;
   }
   //
   for(std::size_t n=0; n < dim; n++) {
      for(std::size_t m=0; m < dim; m++) {
        for(std::size_t i=0; i < size; i++) {
          precompute.push_back( std::real( F[n][i]* F_conj[m][i] ) );
        }
      }
    }
  }

  void update(af::const_ref<FloatType> const& x) {
    MMTBX_ASSERT(x.size() == dim);
    target_ = 0;
    gradient_.fill(0);
    af::shared<FloatType> i_model(i_obs.size());
    i_model.fill(0);
    std::size_t index = 0;
    for(std::size_t n=0; n < dim; n++) {
      for(std::size_t m=0; m < dim; m++) {
        FloatType xnm = x[n]*x[m];
        for(std::size_t i=0; i < size; i++) {
          i_model[i] += xnm * precompute[index];
          index+=1;
        }
      }
    }
    // target
    for(std::size_t i=0; i < size; i++) {
      FloatType diff = i_model[i] - i_obs[i];
      target_ += (diff*diff);
    }
    target_ /= (4*i_obs_sum);
    // grads
    af::shared<FloatType> tmp(i_obs.size());
    index = 0;
    for(std::size_t n=0; n < dim; n++) {
      tmp.fill(0);
      for(std::size_t m=0; m < dim; m++) {
        for(std::size_t i=0; i < size; i++) {
          tmp[i] += x[m] * precompute[index];
          index+=1;
        }
      }
      for(std::size_t i=0; i < size; i++) {
        gradient_[n] += (tmp[i] * ( i_model[i] - i_obs[i] ))/i_obs_sum;
      }
    }
  }

  FloatType target()               { return target_; }
  af::shared<FloatType> gradient() { return gradient_; }

private:
  FloatType target_;
  af::shared<FloatType> gradient_;
  af::shared<af::shared<ComplexType> > F_conj;
  std::size_t dim, size;
  af::shared<af::const_ref<ComplexType> > F;
  af::shared<FloatType> i_obs;
  FloatType i_obs_sum;
  af::shared<FloatType> precompute;
};

/*
Normal matrix A(n,m) = sum_i Re( G_n[i] * conj(G_m[i]) ) for a list of complex
arrays G_n of equal size.
*/
inline
af::versa<double, af::c_grid<2> >
normal_matrix(boost::python::list const& G_)
{
  std::size_t dim = boost::python::len(G_);
  af::shared<af::const_ref<ComplexType> > G(dim);
  for(std::size_t n=0; n < dim; n++) {
    boost::python::extract<af::const_ref<ComplexType> > proxy(G_[n]);
    G[n] = proxy();
    MMTBX_ASSERT(G[n].size() == G[0].size());
  }
  af::versa<double, af::c_grid<2> > A(af::c_grid<2>(dim, dim), 0);
  std::size_t size = (dim == 0) ? 0 : G[0].size();
  for(std::size_t n=0; n < dim; n++) {
    for(std::size_t m=n; m < dim; m++) {
      double s = 0;
      const ComplexType* gn = G[n].begin();
      const ComplexType* gm = G[m].begin();
      for(std::size_t i=0; i < size; i++) {
        s += gn[i].real()*gm[i].real() + gn[i].imag()*gm[i].imag();
      }
      A(n,m) = s;
      A(m,n) = s;
    }
  }
  return A;
}

/*
Solver of the linear system of alg4 with bounds on the unknowns:

  minimize x^T A x / 2 - b^T x  subject to  lower <= x <= upper,

A symmetric positive semi-definite (a normal matrix). Projected Gauss-Seidel:
one unknown at a time is set to its optimal value given the others and
clipped to its bounds. Without active bounds the result is the solution of
A x = b. An unknown with zero diagonal element (empty column) is not
determined by the data; it keeps its starting value.

x: starting values, used as given after clipping to the bounds.
*/
inline
af::shared<double>
alg4_solve_bounded(
  af::const_ref<double, af::c_grid<2> > const& a,
  af::const_ref<double> const& b,
  af::const_ref<double> const& lower,
  af::const_ref<double> const& upper,
  af::const_ref<double> const& x_start,
  int    const& max_sweeps = 10000,
  double const& eps = 1.e-10)
{
  std::size_t dim = b.size();
  MMTBX_ASSERT(a.accessor()[0] == dim && a.accessor()[1] == dim);
  MMTBX_ASSERT(lower.size() == dim && upper.size() == dim);
  MMTBX_ASSERT(x_start.size() == dim);
  af::shared<double> x(dim);
  for(std::size_t i=0; i < dim; i++) {
    MMTBX_ASSERT(lower[i] <= upper[i]);
    x[i] = std::min(upper[i], std::max(lower[i], x_start[i]));
  }
  for(int sweep=0; sweep < max_sweeps; sweep++) {
    double max_shift = 0;
    for(std::size_t i=0; i < dim; i++) {
      double aii = a(i,i);
      double xi = x[i];
      if(aii > 0) {
        double s = b[i];
        for(std::size_t j=0; j < dim; j++) {
          if(j != i) s -= a(i,j)*x[j];
        }
        xi = s/aii;
      }
      xi = std::min(upper[i], std::max(lower[i], xi));
      max_shift = std::max(max_shift, std::abs(xi - x[i]));
      x[i] = xi;
    }
    if(max_shift <= eps) break;
  }
  return x;
}

/*
Phased simultaneous search (alg4) for the coefficients x_n of

  F_model = base + sum_n x_n G_n

that minimize sum_i | F_obs[i] exp(i phi[i]) - F_model[i] |^2, where phi are
the phases of the current F_model, with bounds lower <= x <= upper. Phases
are updated n_cycles times; the normal matrix does not depend on phases and
is computed once.
*/
inline
af::shared<double>
alg4(
  boost::python::list const& G_,
  af::const_ref<double> const& f_obs,
  af::const_ref<ComplexType> const& base,
  af::const_ref<double> const& lower,
  af::const_ref<double> const& upper,
  af::const_ref<double> const& x_start,
  int const& n_cycles = 5,
  af::const_ref<double, af::c_grid<2> > const& A_given =
    af::const_ref<double, af::c_grid<2> >(0, af::c_grid<2>(0,0)))
{
  std::size_t dim = boost::python::len(G_);
  std::size_t size = f_obs.size();
  MMTBX_ASSERT(base.size() == size);
  af::shared<af::const_ref<ComplexType> > G(dim);
  for(std::size_t n=0; n < dim; n++) {
    boost::python::extract<af::const_ref<ComplexType> > proxy(G_[n]);
    G[n] = proxy();
    MMTBX_ASSERT(G[n].size() == size);
  }
  // The normal matrix does not depend on the phases: it can be computed
  // once (normal_matrix) and given here.
  af::versa<double, af::c_grid<2> > A_own;
  af::const_ref<double, af::c_grid<2> > A = A_given;
  if(A.size() == 0 && dim > 0) {
    A_own = normal_matrix(G_);
    A = A_own.const_ref();
  }
  MMTBX_ASSERT(dim == 0 ||
    (A.accessor()[0] == dim && A.accessor()[1] == dim));
  af::shared<double> x(dim);
  for(std::size_t n=0; n < dim; n++) {
    x[n] = std::min(upper[n], std::max(lower[n], x_start[n]));
  }
  af::shared<ComplexType> delta(size);
  af::shared<double> b(dim);
  for(int cycle=0; cycle < n_cycles; cycle++) {
    // F_obs with the phases of the current model, minus the fixed part
    for(std::size_t i=0; i < size; i++) {
      ComplexType m = base[i];
      for(std::size_t n=0; n < dim; n++) m += x[n]*G[n][i];
      double am = std::abs(m);
      ComplexType u = (am > 0) ? m/am : ComplexType(1, 0);
      delta[i] = f_obs[i]*u - base[i];
    }
    for(std::size_t n=0; n < dim; n++) {
      double s = 0;
      const ComplexType* gn = G[n].begin();
      for(std::size_t i=0; i < size; i++) {
        s += gn[i].real()*delta[i].real() + gn[i].imag()*delta[i].imag();
      }
      b[n] = s;
    }
    x = alg4_solve_bounded(A, b.const_ref(), lower, upper,
      x.const_ref());
  }
  return x;
}

}} // namespace mmtbx::mosaic

#endif // MMTBX_MOSAIC_H
