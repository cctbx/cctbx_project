#include <iostream>
#include <cctbx/sgtbx/symbols.h>
#include <cctbx/sgtbx/groups.h>
#include <cctbx/uctbx.h>

int main()
{
  uctbx::UnitCell UnitCell(uctbx::uc_params(11, 12, 13, 90, 100, 90));
  std::cout << UnitCell << std::endl;
  sgtbx::SpaceGroupSymbols Symbols("C 2");
  sgtbx::SpaceGroup SgOps(Symbols.Hall());
  SgOps.CheckUnitCell(UnitCell);
  for (int i = 0; i < SgOps.OrderZ(); i++) {
    std::cout << SgOps(i).as_xyz() << std::endl;
  }
  return 0;
}
