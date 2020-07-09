#ifndef MMTBX_MOSAIC_H
#define MMTBX_MOSAIC_H

#include <cfloat>
#include <mmtbx/error.h>
#include <mmtbx/import_scitbx_af.h>
#include <cctbx/sgtbx/space_group.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <mmtbx/f_model/f_model.h>
#include <scitbx/matrix/outer_product.h>
#include <scitbx/array_family/versa_algebra.h>
#include <scitbx/array_family/shared_algebra.h>
#include <scitbx/array_family/small_algebra.h>
#include <scitbx/matrix/eigensystem.h>
#include <scitbx/matrix/packed.h>
#include <scitbx/array_family/versa_matrix.h>
#include <scitbx/math/cubic_equation.h>
#include <mmtbx/bulk_solvent/bulk_solvent.h>
#include <boost/python/list.hpp>

namespace mmtbx { namespace mosaic {


template <typename FloatType, typename ComplexType>
 af::shared<FloatType>
 alg4(
   boost::python::list        const& F_,
   af::const_ref<ComplexType> const& f_obs,
   FloatType                  const& k)
 {
   int dim = boost::python::len(F_);
   af::shared<af::shared<FloatType> > F(dim);
   for(std::size_t i=0; i < dim; i++) {
     boost::python::extract<af::shared<double> > elem_proxy_1(F_[i]);
     F.push_back(elem_proxy_1());
   }
   int size = F[0].size();
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
     }
     FloatType Hj = 0;
     for(std::size_t i=0; i < size; i++) {
       Hj += std::real(Fj[i]*std::conj(f_obs[i]));
     }
     b[j] = Hj;
   }
   af::versa<FloatType, af::c_grid<2> > A_inv(
      scitbx::matrix::packed_u_as_symmetric(
        scitbx::matrix::eigensystem::real_symmetric<FloatType>(
          A.const_ref(), /*relative_epsilon*/ 1.e-9,/*absolute_epsilon*/ 1.e-9)
            .generalized_inverse_as_packed_u().const_ref()
            )
            );
   af::shared<FloatType> x = af::matrix_multiply(
     A_inv.const_ref(), b.const_ref());
   return x;


   // LSE example
   //std::size_t n = 2;
   //af::versa<FloatType, af::mat_grid> m(af::mat_grid(n, n), 0);
   //m(0,0)=1;
   //m(0,1)=3;
   //m(1,0)=3;
   //m(1,1)=-1;
   //af::shared<FloatType> b(n);
   //b[0] = 5;
   //b[1] = 3;
   //
   //af::versa<FloatType, af::c_grid<2> > m_inv(
   //   scitbx::matrix::packed_u_as_symmetric(
   //     scitbx::matrix::eigensystem::real_symmetric<FloatType>(
   //       m.const_ref(), /*relative_epsilon*/ 1.e-9,/*absolute_epsilon*/ 1.e-9)
   //         .generalized_inverse_as_packed_u().const_ref()
   //         )
   //         );
   //
   //af::shared<FloatType> x = af::matrix_multiply(
   //  m_inv.const_ref(), b.const_ref());
   //
   //std::cout<<x[0]<<" "<<x[1]<< std::endl;
   //
   //return b;
 };


}} // namespace mmtbx::mosaic

#endif // MMTBX_MOSAIC_H
