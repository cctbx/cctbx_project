#include <iostream>
#include <cctbx/sgtbx/symbols.h>
#include <cctbx/sgtbx/groups.h>
#include <cctbx/uctbx.h>

int main()
{
  cctbx::uctbx::UnitCell UnitCell(cctbx::uctbx::uc_params(11,12,13,90,100,90));
  std::cout << UnitCell << std::endl;
  cctbx::sgtbx::SpaceGroupSymbols Symbols("C 2");
  cctbx::sgtbx::SpaceGroup SgOps(Symbols);
  SgOps.CheckUnitCell(UnitCell);
  for (int i = 0; i < SgOps.OrderZ(); i++) {
    std::cout << SgOps(i).as_xyz() << std::endl;
  }
  return 0;
}
