// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <iostream>
#include <cctbx/sgtbx/groups.h>
#include <cctbx/array_family/simple_io.h>

int main(void)
{
  using namespace cctbx::sgtbx;
  parse_string HSym("P 31");
  SpaceGroup SgOps;
  try {
    SgOps.ParseHallSymbol(HSym);
  }
  catch (const cctbx::error& e) {
    std::cout << e.what() << std::endl;
    std::cout << HSym.string() << std::endl;
    for (int i = 0; i < HSym.where(); i++) std::cout << "_";
    std::cout << "^" << std::endl;
  }

  RTMx M = SgOps(1);
  std::cout << M.as_array(double(0)).ref() << std::endl;
  std::cout << M.as_int_array().ref() << std::endl;

  return 0;
}
