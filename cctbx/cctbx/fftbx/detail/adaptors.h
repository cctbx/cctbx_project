// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Dec 21: iterator-based interface (rwgk)
     2001 Nov 03: fftbx started, based on fftpack41 (rwgk)
 */

#ifndef CCTBX_FFTBX_DETAIL_ADAPTORS_H
#define CCTBX_FFTBX_DETAIL_ADAPTORS_H

#include <boost/array.hpp>
#include <cctbx/ndim.h>

namespace cctbx { namespace fftbx {
  namespace detail {

    template <std::size_t D,
              typename IteratorType,
              typename ValueType = typename IteratorType::value_type>
    class access_tp
    {
      public:
        typedef IteratorType iterator_type;
        typedef ValueType value_type;

        access_tp(iterator_type Start,
                  const std::size_t& N0,
                  const std::size_t& N1)
          : m_Start(Start) {
          m_N[0] = N0;
          m_N[1] = N1;
        }
        access_tp(iterator_type Start,
                  const std::size_t& N0,
                  const std::size_t& N1,
                  const std::size_t& N2)
          : m_Start(Start) {
          m_N[0] = N0;
          m_N[1] = N1;
          m_N[2] = N2;
        }
        value_type&
        operator()(const std::size_t& i0,
                   const std::size_t& i1) {
          boost::array<std::size_t, D> I;
          I[0] = i0;
          I[1] = i1;
          return m_Start[fortran_index_1d<D>()(m_N, I)];
        }
        value_type&
        operator()(const std::size_t& i0,
                   const std::size_t& i1,
                   const std::size_t& i2) {
          boost::array<std::size_t, D> I;
          I[0] = i0;
          I[1] = i1;
          I[2] = i2;
          return m_Start[fortran_index_1d<D>()(m_N, I)];
        }
      private:
        iterator_type m_Start;
        boost::array<std::size_t, D> m_N;
    };

  } // namespace detail
}} // namespace cctbx::fftbx

#endif // CCTBX_FFTBX_DETAIL_ADAPTORS_H
