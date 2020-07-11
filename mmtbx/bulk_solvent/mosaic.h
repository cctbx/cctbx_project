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

// Phased simultaneous search (alg4). Python equivalent exists.
template <typename FloatType>
 af::shared<FloatType>
 alg4(
   boost::python::list      const& F_,
   af::const_ref<FloatType> const& f_obs,
   int                      const& max_cycles,
   FloatType                const& auto_converge_eps)
 {
   int dim = boost::python::len(F_);
   // Extract Fs
   af::shared<af::const_ref<ComplexType> > F(dim);
   for(std::size_t i=0; i < dim; i++) {
     boost::python::extract<af::const_ref<ComplexType> > elem_proxy_1(F_[i]);
     af::const_ref<ComplexType>  fm = elem_proxy_1();
     F[i] = fm;
   }
   //
   int size = f_obs.size();
   // Get copy of F[0] which Fc
   af::shared<ComplexType> fc_d(F[0].size());
   for(std::size_t i=0; i < size; i++) fc_d[i]=F[0][i];
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
   af::shared<FloatType> x_res(dim, 0);
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
       af::const_ref<ComplexType> Fj;
       if(j==0) Fj = fc_d.const_ref();
       else     Fj = F[j];
       for(std::size_t n=0; n < dim; n++) {
         FloatType Gjn = 0;
         for(std::size_t i=0; i < size; i++) {
           if(n==0) Gjn += std::real(Fj[i]*std::conj(fc_d[i]) );
           else     Gjn += std::real(Fj[i]*F_CONJ[n][i] );
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
     for(std::size_t i=0; i < dim; i++) {
       if(i==0) continue;
       for(std::size_t j=0; j < size; j++) {
         fc_d[j] += x[i] * F[i][j];
       }
     }
     x_res = x_res + x;
     x_[0]=x[0];
     for(std::size_t i=1; i < dim; i++) x_[i] = x_res[i];
     //
     cntr+=1;
     // Exit endless loop conditions
     if(cntr>max_cycles) break;
     if(cntr==0) { for(std::size_t i=0; i < dim; i++) x_prev[i] = x_[i]; }
     else {
       FloatType max_diff = std::abs(x_prev[0]-x_[0]);
       for(std::size_t i=1; i < dim; i++) {
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
