// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <cmath>
#include <cctbx/sgtbx/groups.h>
#include <cctbx/sgtbx/coordinates.h>
#include <cctbx/basic/define_range.h>

namespace sgtbx {

  namespace {

    uctbx::Vec3 ModShort(uctbx::Vec3 X)
    {
      uctbx::Vec3 result;
      rangei(3) {
        result[i] = std::fmod(X[i], 1.);
        if (result[i] < -.5)
          result[i] += 1.;
        else if (result[i] > .5)
          result[i] -= 1.;
      }
      return result;
    }

    double getDelta2(const uctbx::UnitCell& uc,
                     const uctbx::Vec3& X, const uctbx::Vec3& Y)
    {
      uctbx::Vec3 Delta = ModShort(X - Y);
      uctbx::Vec3 CartDelta = uc.orthogonalize(Delta);
      return CartDelta * CartDelta;
    }
  }

  double
  SymEquivCoordinates::getMinDelta2(const uctbx::UnitCell& uc,
                                    const uctbx::Vec3& X) const
  {
    double MinDelta2 = getDelta2(uc, X, m_Coordinates[0]);
    for (int i = 1; i < m_Coordinates.size(); i++) {
      double Delta2 = getDelta2(uc, X, m_Coordinates[i]);
      if (MinDelta2 > Delta2)
          MinDelta2 = Delta2;
    }
    return MinDelta2;
  }

  SymEquivCoordinates::SymEquivCoordinates(const uctbx::UnitCell& uc,
                                           const SgOps& sgo,
                                           const uctbx::Vec3& X,
                                           double MinimumDistance)
    : m_Coordinates(1, X)
  {
    double MD2 = MinimumDistance * MinimumDistance;
    rangei(sgo.OrderZ()) {
      uctbx::Vec3 SX = sgo(i) * X;
      if (getMinDelta2(uc, SX) >= MD2)
        m_Coordinates.push_back(SX);
    }
    cctbx_assert(sgo.OrderZ() % m_Coordinates.size() == 0);
  }

  std::complex<double>
  SymEquivCoordinates::StructureFactor(const Miller::Index& H) const
  {
    using cctbx::constants::pi;
    std::complex<double> F(0., 0.);
    rangei(M()) {
      double phase = 2. * pi * (H * m_Coordinates[i]);
      F += std::complex<double>(std::cos(phase), std::sin(phase));
    }
    return F;
  }

} // namespace sgtbx
