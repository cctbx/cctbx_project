// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <iostream>
#include <cctbx/sgtbx/symbols.h>
#include <cctbx/sgtbx/groups.h>
#include <cctbx/uctbx.h>

int main()
{
  uctbx::UnitCell UnitCell(uctbx::uc_params(11, 12, 13, 90, 100, 90));
  std::cout << UnitCell << std::endl;
  sgtbx::SpaceGroupSymbols Symbols("C 2");
  sgtbx::SgOps SgOps(Symbols.Hall());
  SgOps.CheckUnitCell(UnitCell);
  for (int i = 0; i < SgOps.OrderZ(); i++) {
    std::cout << SgOps(i).as_xyz() << std::endl;
  }
  return 0;
}
