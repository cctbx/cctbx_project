// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_SGTBX_COORDINATES_H
#define CCTBX_SGTBX_COORDINATES_H

#include <vector>
#include <complex>
#include <cctbx/uctbx.h>
#include <cctbx/miller.h>
#include <cctbx/sgtbx/matrix.h>

namespace sgtbx {

  inline uctbx::Vec3 operator*(const RTMx& lhs, const uctbx::Vec3& rhs) {
    const RotMx& R = lhs.Rpart();
    const TrVec& T = lhs.Tpart();
    double TBF = T.BF();
    uctbx::Vec3 result;
    result[0] = R[0] * rhs[0] + R[1] * rhs[1] + R[2] * rhs[2] + T[0] / TBF;
    result[1] = R[3] * rhs[0] + R[4] * rhs[1] + R[5] * rhs[2] + T[1] / TBF;
    result[2] = R[6] * rhs[0] + R[7] * rhs[1] + R[8] * rhs[2] + T[2] / TBF;
    return result;
  }

  inline double operator*(const uctbx::Vec3& lhs, const uctbx::Vec3& rhs) {
    double result = 0.;
    for(int i=0;i<3;i++) result += lhs[i] * rhs[i];
    return result;
  }

  inline uctbx::Vec3 operator-(const uctbx::Vec3& lhs, const uctbx::Vec3& rhs)
  {
    uctbx::Vec3 result;
    for(int i=0;i<3;i++) result[i] = lhs[i] - rhs[i];
    return result;
  }

  inline double operator*(const Miller::Index& lhs, const uctbx::Vec3& rhs) {
    double result = 0.;
    for(int i=0;i<3;i++) result += lhs[i] * rhs[i];
    return result;
  }

  //! Container for symmetry equivalent (atomic) coordinates.
  class SymEquivCoordinates {
    public:
      //! Compute list of symmetry equivalent coordinates.
      /*! The symmetry operations (SgOps) are applied to the
          fractional coordinates X. The unit cell parameters (UnitCell)
          are used to compute the distances between the symmetry
          equivalent coordinates. If the distance between a symmetry
          equivalent coordinate and the input coordinate X is greater
          than or equal to the MinimumDistance, the symmetry equivalent
          coordinate is added to the list. Otherwise X is a special
          position, and the symmetry equivalent coordinate is not added
          to the list.
          <p>
          The MinimumDistance should be large enough to account for
          rounding errors. However, the tolerance should not be too
          large because the current implementation does not move
          the X coordinate to the exact special position.
          <p>
          The simple distance calculation used in the current
          implementation is not guaranteed to be numerically stable.
          However, in general this should not be a problem. As a
          safeguard it is asserted that the number of symmetry
          equivalent coordinates (M()) in the list is a factor of the
          space group multiplicity.
       */
      SymEquivCoordinates(const uctbx::UnitCell& uc,
                          const SgOps& sgo,
                          const uctbx::Vec3& X,
                          double MinimumDistance = 0.05);
      //! Number of symmetry equivalent coordinates (multiplicity).
      inline int M() const { return m_Coordinates.size(); }
      //! Return the i'th symmetry equivalent coordinate.
      /*! An exception is thrown if i is out of range.
       */
      inline const uctbx::Vec3& operator()(int i) const {
        if (i < 0 || i > M()) throw error("Index out of range.");
        return m_Coordinates[i];
      }
      //! Compute Sum(exp(2 pi i H X)) for all symmetry equivalent X.
      /*! This sum is a sub-expression in the structure factor
          calculation. See file examples/python/generate_hklf.py.
       */
      std::complex<double> StructureFactor(const Miller::Index& H) const;

    private:
      std::vector<uctbx::Vec3> m_Coordinates;

      double
      getMinDelta2(const uctbx::UnitCell& uc, const uctbx::Vec3& X) const;
  };

} // namespace sgtbx

#endif // CCTBX_SGTBX_COORDINATES_H
