// $Id$

#include <iostream>
#include <cctbx/sgtbx/symbols.h>
#include <cctbx/sgtbx/groups.h>
#include <cctbx/uctbx.h>

int main()
{
  uctbx::UnitCell UnitCell(uctbx::uc_params(10, 10, 15, 90, 90, 120));
  std::cout << UnitCell << std::endl;
  sgtbx::SpaceGroupSymbols Symbols("P 62 2 2");
  sgtbx::SgOps SgOps(Symbols.Hall());
  SgOps.CheckUnitCell(UnitCell);
  for (int i = 0; i < SgOps.OrderZ(); i++) {
    std::cout << SgOps(i).as_xyz() << std::endl;
  }
  return 0;
}
