// $Id$

#include <iostream>
#include <cctbx/sgtbx/groups.h>

int main(void)
{
  using namespace sgtbx;
  parse_string HSym("P x");
  SgOps sgops;
  try {
    sgops.ParseHallSymbol(HSym);
  }
  catch (const error& e) {
    std::cout << e.what() << std::endl;
    std::cout << HSym.string() << std::endl;
    for (int i = 0; i < HSym.where(); i++) std::cout << "_";
    std::cout << "^" << std::endl;
  }

  return 0;
}
