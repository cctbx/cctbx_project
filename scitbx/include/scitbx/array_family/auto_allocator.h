/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copied from cctbx/array_family (R.W. Grosse-Kunstleve)
     2002 Feb: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_ARRAY_FAMILY_AUTO_ALLOCATOR_H
#define SCITBX_ARRAY_FAMILY_AUTO_ALLOCATOR_H

#include <cstddef>
#include <boost/config.hpp>

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

#include <scitbx/array_family/alignment_calculator.h>

namespace scitbx { namespace af { namespace detail {

  template <typename T, std::size_t N>
  union auto_allocator {
    char buffer[N * sizeof(T)];
    typename scitbx::af::alignment_calculator::AlignmentCalculator<
      SCITBX_AF_TYPELIST1(T)>::Result dummy_;
  };

}}} // namespace scitbx::af::detail

#else // !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

namespace scitbx { namespace af { namespace detail {

  template <typename T>
  struct SimpleAlignmentCalculator { typedef long double Result; };

#define SCITBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(T) \
template <> struct SimpleAlignmentCalculator<T> { typedef T Result; };

  SCITBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(char)
  SCITBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(signed char)
  SCITBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(unsigned char)
  SCITBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(short int)
  SCITBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(unsigned short int)
  SCITBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(int)
  SCITBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(unsigned int)
  SCITBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(long int)
  SCITBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(unsigned long int)
  SCITBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(float)
  SCITBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(double)
  SCITBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(long double)
  SCITBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(char*)
  SCITBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(short int*)
  SCITBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(int*)
  SCITBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(long int*)
  SCITBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(float*)
  SCITBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(double*)
  SCITBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(long double*)
  SCITBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR(void*)

#undef SCITBX_ARRAY_FAMILY_DEFINE_SIMPLE_ALIGNMENT_CALCULATOR

  template <typename T, std::size_t N>
  union auto_allocator {
    char buffer[N * sizeof(T)];
    typename SimpleAlignmentCalculator<T>::Result dummy_;
  };

}}} // namespace scitbx::af::detail

#endif // !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

#endif // SCITBX_ARRAY_FAMILY_AUTO_ALLOCATOR_H
