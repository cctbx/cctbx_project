/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Dec 2001: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_VECTOR_UNARY_H
#define CCTBX_VECTOR_UNARY_H

#include <complex>
#include <cctbx/shared_storage.h>
#include <cctbx/basic/meta.h>

namespace cctbx { namespace vector {

  template <typename UnaryOperation,
            typename ArgumentVectorType,
            typename ResultValueType>
  shared_storage<ResultValueType>
  apply(UnaryOperation op,
        const ArgumentVectorType& arguments,
        type_holder<ResultValueType> result_type_holder)
  {
    shared_storage<ResultValueType> result(arguments.size());
    for (std::size_t i=0;i<arguments.size();i++) {
      result[i] = op(arguments[i]);
    }
    return result;
  }

  template <typename UnaryOperation,
            typename ArgumentVectorType>
  shared_storage<typename UnaryOperation::result_type>
  apply(UnaryOperation op,
        const ArgumentVectorType& arguments)
  {
    return apply(op, arguments,
                 type_holder<typename UnaryOperation::result_type>());
  }

  template <typename ComplexVectorType, typename ResultVectorType>
  ResultVectorType complex_abs(const ComplexVectorType& v,
                               type_holder<ResultVectorType>)
  {
    ResultVectorType result(v.size());
    for(std::size_t i=0;i<v.size();i++)
      result[i] = v[i].real() * v[i].real() + v[i].imag() * v[i].imag();
    return result;
  }

}} // namespace cctbx::vector

#endif // CCTBX_VECTOR_UNARY_H
