#ifndef CCTBX_SGTBX_ASU_SHORTCUTS_H
#define CCTBX_SGTBX_ASU_SHORTCUTS_H

#include "cut.h"

namespace cctbx { namespace sgtbx { namespace asu {

  // shortcuts
  const rational_t r1 = 1;
  const cut
    x1 = cut(-1,0,0, r1),
    x0 = -x1*0,
    x2 = x1/2,
    x3 = x1/3,
    x4 = x1/4,
    x8 = x1/8,
    y1 = cut(0,-1,0, r1),
    y0 = -y1*0,
    y2 = y1/2,
    y3 = y1/3,
    y4 = y1/4,
    y8 = y1/8,
    z1 = cut(0,0,-1, r1),
    z0 = -z1*0,
    z2 = z1/2,
    z3 = z1/3,
    z4 = z1/4,
    z6 = z1/6,
    z8 = z1/8,
    z12 = z1/12,
    p1 = cut(-1,1,0, r1),
    p0 = -p1*0,
    p2 = p1/2,
    p3 = p1/3,
    p4 = p1/4,
    m1 = cut(-1,-1,0, r1),
    m0 = -m1*0,
    m2 = m1/2,
    h1 = cut(1,-2,0, r1),
    h0 = -h1*0,
    k1 = cut(-2,1,0, r1),
    k0 = -k1*0
  ;

}}}
#endif
