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

#include <iostream>
#include <cctbx/fixes/cmath>
#include <cctbx/fixes/cstdlib>
#include <boost/array.hpp>
#include <cctbx/basic/matrixlite.h>
#include <cctbx/constants.h>

namespace cctbx {
  //! Miller index namespace.
  namespace Miller {

    typedef MatrixLite::itype::Vec3 Vec3;

    //! Enumeration for symbolic subscripting (e.g. MillerIndex[H]).
    enum {H, K, L};

    //! Miller index class.
    class Index : public Vec3 {
      public:
        //! @name Constructors.
        //@{
        Index() { for(int i=0;i<3;i++) elems[i] = 0; }
        Index(const int *hkl) { for(int i=0;i<3;i++) elems[i] = hkl[i]; }
        Index(const int h, const int k, const int l) {
          elems[0] = h; elems[1] = k; elems[2] = l;
        }
        //@}

        //! @name Convenience methods.
        //@{
        inline bool is000() const {
          return !(elems[0] || elems[1] || elems[2]);
        }
        inline Index operator-() const {
          return Index(-elems[0], -elems[1], -elems[2]);
        }
        inline Index FriedelMate() const {
          return operator-();
        }
        //@}

        //! @name Test for equality and inequality.
        //@{
        inline friend bool operator==(const Index& lhs, const Index& rhs) {
          for(int i=0;i<3;i++) if (lhs[i] != rhs[i]) return false;
          return true;
        }
        inline friend bool operator!=(const Index& lhs, const Index& rhs) {
          return !(lhs == rhs);
        }
        //@}

        //! @name Definition of sort order for human-readable listings.
        //@{
        /*! This comparison is computationally more expensive than
            the Miller::hashCompare below.
         */
        inline bool operator<(const Index& m2) const
        {
          const int P[3] = {2, 0, 1};
          int i;
          for(i=0;i<3;i++) {
            if (elems[P[i]] >= 0 && m2[P[i]] <  0) return true;
            if (elems[P[i]] <  0 && m2[P[i]] >= 0) return false;
          }
          for(i=0;i<3;i++) {
            if (std::abs(elems[P[i]]) < std::abs(m2[P[i]])) return true;
            if (std::abs(elems[P[i]]) > std::abs(m2[P[i]])) return false;
          }
          return false;
        }
        inline bool operator>(const Index& m2) const {
          return !(*this < m2);
        }
        //@}
    };

    //! Multiplication of Miller indices and fractional coordiantes.
    inline double
    operator*(const Index& lhs, const MatrixLite::dtype::Vec3& rhs) {
      double result = 0.;
      for(int i=0;i<3;i++) result += lhs[i] * rhs[i];
      return result;
    }

    //! Definition of fast comparison for use in, e.g., std::map<Miller::Index>.
    class hashCompare {
      public:
        //! This fast comparison function is implemented as operator().
        inline bool operator()(const Index& m1,const Index& m2) const {
          for(int i=0;i<3;i++) {
            if (m1[i] < m2[i]) return true;
            if (m1[i] > m2[i]) return false;
          }
          return false;
        }
    };

    //! iostream output operator for class Miller::Index.
    inline std::ostream& operator<<(std::ostream& os, const Miller::Index& MIx)
    {
      os << "H=" << MIx[0] << " K=" << MIx[1] << " L=" << MIx[2];
      return os;
    }

  } // namespace Miller
} // namespace cctbx

#endif // CCTBX_MILLER_H
