#ifndef MMTBX_MOSAIC_H
#define MMTBX_MOSAIC_H

#include <mmtbx/import_scitbx_af.h>
#include <scitbx/matrix/eigensystem.h>
#include <scitbx/array_family/versa_matrix.h>
#include <boost/python/list.hpp>

namespace mmtbx { namespace mosaic {


template <typename FloatType, typename ComplexType>
 af::shared<FloatType>
 alg4(
   boost::python::list        const& F_,
   af::const_ref<ComplexType> const& f_obs)
 {
   int dim = boost::python::len(F_);
   af::shared<af::shared<ComplexType> > F(dim);
   for(std::size_t i=0; i < dim; i++) {
     boost::python::extract<af::shared<ComplexType> > elem_proxy_1(F_[i]);
     af::shared<ComplexType> fm = elem_proxy_1();
     F[i] = fm;
   }
   int size = f_obs.size();

   //af::shared<ComplexType> fc_d = F[0];
   //af::shared<FloatType> x_res(dim, 0);
   //af::shared<FloatType> x_prev(dim, 0);
   //int cntr = 0;
   //while(true) {

     af::shared<FloatType> b(dim);
     af::versa<FloatType, af::mat_grid> A(af::mat_grid(dim, dim), 0);
     for(std::size_t j=0; j < dim; j++) {
       af::shared<ComplexType> Fj = F[j];
       for(std::size_t n=0; n < dim; n++) {
         af::shared<ComplexType> Fn = F[n];
         FloatType Gjn = 0;
         for(std::size_t i=0; i < size; i++) {
           Gjn += std::real(Fj[i]*std::conj(Fn[i]));
         }
         A(j,n) = Gjn;
       }
       FloatType Hj = 0;
       for(std::size_t i=0; i < size; i++) {
         Hj += std::real(Fj[i]*std::conj(f_obs[i]));
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
     //
     
   //  for(std::size_t i=0; i < dim; i++) {
   //    if(i==0) continue;
   //    x[i]
   //  }
   //
   //}

   return x;
 };


}} // namespace mmtbx::mosaic

#endif // MMTBX_MOSAIC_H
