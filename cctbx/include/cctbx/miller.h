// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Jul 02: Merged from CVS branch sgtbx_special_pos (rwgk)
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

/*! \file
    Toolbox for Miller indices.
 */

#ifndef CCTBX_MILLER_H
#define CCTBX_MILLER_H

#if defined(__GNUC__) && __GNUC__ < 3
# include <iostream>
#else
# include <ostream>
#endif
#include <cctbx/fixes/cstdlib>
#include <cctbx/error.h>
#include <cctbx/coordinates.h>
#include <cctbx/array_family/tiny_types.h>
#include <cctbx/math/utils.h>

namespace cctbx {
  //! Miller index namespace.
  namespace miller {

    //! Enumeration for symbolic subscripting (e.g. MillerIndex[H]).
    enum {H, K, L};

    //! Miller index class.
    class Index : public af::int3 {
      public:
        //! @name Constructors.
        //@{
        Index() {
          for(std::size_t i=0;i<3;i++) elems[i] = 0;
        }
        Index(const af::int3& v) {
          for(std::size_t i=0;i<3;i++) elems[i] = v[i];
        }
        explicit Index(const int* hkl) {
          for(std::size_t i=0;i<3;i++) elems[i] = hkl[i];
        }
        Index(int h, int k, int l) {
          elems[0] = h; elems[1] = k; elems[2] = l;
        }
        //@}

        //! @name Convenience methods.
        //@{
        bool is000() const {
          return !(elems[0] || elems[1] || elems[2]);
        }
        Index operator-() const {
          return Index(-elems[0], -elems[1], -elems[2]);
        }
        //@}

        //! @name Definition of sort order for human-readable listings.
        //@{
        /*! This comparison is computationally more expensive than
            the miller::hashCompare below.
         */
        bool operator<(const Index& m2) const
        {
          const int P[3] = {2, 0, 1};
          std::size_t i;
          for(i=0;i<3;i++) {
            if (elems[P[i]] >= 0 && m2[P[i]] <  0) return true;
            if (elems[P[i]] <  0 && m2[P[i]] >= 0) return false;
          }
          for(i=0;i<3;i++) {
            if (math::abs(elems[P[i]]) < math::abs(m2[P[i]])) return true;
            if (math::abs(elems[P[i]]) > math::abs(m2[P[i]])) return false;
          }
          return false;
        }
        bool operator>(const Index& m2) const
        {
          if (*this < m2) return false;
          for(std::size_t i=0;i<3;i++) if (elems[i] != m2[i]) return true;
          return false;
        }
        //@}

        //! Test for equality.
        bool operator==(const Index& m2) const
        {
          return af::cmp(this->const_ref(), m2.const_ref()) == 0;
        }

        //! Test for inequality.
        bool operator!=(const Index& m2) const
        {
          return af::cmp(this->const_ref(), m2.const_ref()) != 0;
        }
    };

    //! Multiplication of Miller indices and fractional coordiantes.
    template <class FloatType>
    inline FloatType
    operator*(const Index& lhs, const fractional<FloatType>& rhs) {
      FloatType result = 0.;
      for(std::size_t i=0;i<3;i++) result += lhs[i] * rhs[i];
      return result;
    }

    /*! \brief Definition of fast comparison for use in,
        e.g., std::map<miller::Index>.
     */
    class hashCompare {
      public:
        //! This fast comparison function is implemented as operator().
        bool operator()(const Index& m1,const Index& m2) const {
          for(std::size_t i=0;i<3;i++) {
            if (m1[i] < m2[i]) return true;
            if (m1[i] > m2[i]) return false;
          }
          return false;
        }
    };

    //! ostream output operator for class miller::Index.
    inline std::ostream& operator<<(std::ostream& os, const miller::Index& MIx)
    {
      os << "H=" << MIx[0] << " K=" << MIx[1] << " L=" << MIx[2];
      return os;
    }

  } // namespace miller
} // namespace cctbx

#endif // CCTBX_MILLER_H
