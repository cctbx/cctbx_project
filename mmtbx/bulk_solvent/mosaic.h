#ifndef MMTBX_MOSAIC_H
#define MMTBX_MOSAIC_H

#include <mmtbx/import_scitbx_af.h>
#include <scitbx/matrix/eigensystem.h>
#include <scitbx/array_family/versa_matrix.h>
#include <boost/python/list.hpp>
#include <scitbx/math/utils.h>
#include <mmtbx/error.h>

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

// Phased simultaneous search (alg4). Python equivalent exists.
template <typename FloatType>
 af::shared<FloatType>
 alg4(
   boost::python::list        const& F_,
   af::const_ref<FloatType>   const& f_obs,
   af::const_ref<ComplexType> const& phase_source,
   int                        const& max_cycles,
   FloatType                  const& auto_converge_eps)
 {
   int dim = boost::python::len(F_);
   // Extract Fs
   af::shared<af::const_ref<ComplexType> > F(dim);
   for(std::size_t i=0; i < dim; i++) {
     boost::python::extract<af::const_ref<ComplexType> > elem_proxy_1(F_[i]);
     af::const_ref<ComplexType>  fm = elem_proxy_1();
     F[i] = fm;
     MMTBX_ASSERT(fm.size() == f_obs.size());
     MMTBX_ASSERT(fm.size() == phase_source.size());
   }
   //
   int size = f_obs.size();
   // Get copy of F[0] which Fc
   af::shared<ComplexType> fc_d(phase_source.size());
   for(std::size_t i=0; i < size; i++) fc_d[i]=phase_source[i];
   // Pre-compute complex conj
   af::shared<af::shared<ComplexType> > F_CONJ(dim);
   for(std::size_t i=0; i < dim; i++) {
     af::shared<ComplexType> f(size);
     for(std::size_t j=0; j < size; j++) {
       f[j]=std::conj(F[i][j]);
     }
     F_CONJ[i] = f;
   }
   // Refinement iterations
   af::shared<FloatType> x_prev(dim, 0);
   af::shared<FloatType> x_(dim, 0);
   int cntr = 0;
   while(true) {
     af::shared<FloatType> b(dim);
     af::versa<FloatType, af::mat_grid> A(af::mat_grid(dim, dim), 0);
     // Transfer phases from current F[0] (Fc) onto Fobs
     af::shared<ComplexType> fo_cmpl_conj(size);
     for(std::size_t i=0; i < size; i++) {
       fo_cmpl_conj[i] = std::conj(
         f_obs[i] * scitbx::math::unit_complex(std::arg(fc_d[i])));
     }
     // Construct A and b of the system of linear equations A*x=b
     for(std::size_t j=0; j < dim; j++) {
       af::const_ref<ComplexType> Fj = F[j];
       for(std::size_t n=0; n < dim; n++) {
         FloatType Gjn = 0;
         for(std::size_t i=0; i < size; i++) {
           Gjn += std::real(Fj[i]*F_CONJ[n][i] );
         }
         A(j,n) = Gjn;
       }
       FloatType Hj = 0;
       for(std::size_t i=0; i < size; i++) {
         Hj += std::real(Fj[i]*fo_cmpl_conj[i]);
       }
       b[j] = Hj;
     }
     // Solve SoLE A*x=b: x = A_inv*b
     af::versa<FloatType, af::c_grid<2> > A_inv(
        scitbx::matrix::packed_u_as_symmetric(
          scitbx::matrix::eigensystem::real_symmetric<FloatType>(
            A.const_ref(), /*relative_epsilon*/ 1.e-9,/*absolute_epsilon*/ 1.e-9)
              .generalized_inverse_as_packed_u().const_ref()
              )
              );
     af::shared<FloatType> x = af::matrix_multiply(
       A_inv.const_ref(), b.const_ref());
     // Update F[0] (Fcalc)
     fc_d.fill(0);
     for(std::size_t i=0; i < dim; i++) {
       //if(i==0) continue;
       for(std::size_t j=0; j < size; j++) {
         fc_d[j] += x[i] * F[i][j];
       }
     }
     for(std::size_t i=0; i < dim; i++) x_[i] = x[i];
     //
     cntr+=1;
     // Exit endless loop conditions
     if(cntr>max_cycles) break;
     if(cntr==0) { for(std::size_t i=0; i < dim; i++) x_prev[i] = x_[i]; }
     else {
       FloatType max_diff = std::abs(x_prev[0]-x_[0]);
       for(std::size_t i=0; i < dim; i++) {
         FloatType tmp = std::abs( x_prev[i]-x_[i] );
         if(tmp>max_diff) max_diff=tmp;
       }

       if(max_diff<=auto_converge_eps) break;
       for(std::size_t i=0; i < dim; i++) x_prev[i] = x_[i];
     }
   } // end while loop
   return x_;
 };


}} // namespace mmtbx::mosaic

#endif // MMTBX_MOSAIC_H
