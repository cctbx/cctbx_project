// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Feb 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_AUTO_ALLOCATOR_H
#define CCTBX_ARRAY_FAMILY_AUTO_ALLOCATOR_H

#include <cstddef>
#include <boost/config.hpp>

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

#include <cctbx/array_family/alignment_calculator.h>

namespace cctbx { namespace af { namespace detail {

  template <typename T, std::size_t N>
  union auto_allocator {
    char buffer[N * sizeof(T)];
    typename cctbx::af::alignment_calculator::AlignmentCalculator<
      CCTBX_AF_TYPELIST1(T)>::Result dummy_;
  };

}}} // namespace cctbx::af::detail

#else // !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

namespace cctbx { namespace af { namespace detail {

  template <typename T>
  struct SimpleAlignmentCalculator { typedef long double Result; };

#define CCTBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(T) \
template <> struct SimpleAlignmentCalculator<T> { typedef T Result; };

  CCTBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(char)
  CCTBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(signed char)
  CCTBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(unsigned char)
  CCTBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(short int)
  CCTBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(unsigned short int)
  CCTBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(int)
  CCTBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(unsigned int)
  CCTBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(long int)
  CCTBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(unsigned long int)
  CCTBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(float)
  CCTBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(double)
  CCTBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(long double)
  CCTBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(char*)
  CCTBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(short int*)
  CCTBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(int*)
  CCTBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(long int*)
  CCTBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(float*)
  CCTBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(double*)
  CCTBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(long double*)
  CCTBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(void*)

#undef CCTBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR

  template <typename T, std::size_t N>
  union auto_allocator {
    char buffer[N * sizeof(T)];
    typename SimpleAlignmentCalculator<T>::Result dummy_;
  };

}}} // namespace cctbx::af::detail

#endif // !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

#endif // CCTBX_ARRAY_FAMILY_AUTO_ALLOCATOR_H
