// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <cctbx/uctbx.h>
#include <cctbx/miller.h>

int main(void)
{
  using namespace cctbx;
  uctbx::UnitCell uc(uctbx::uc_params(2, 3, 4, 81, 82, 83));
  std::cout << uc << std::endl;
  uctbx::uc_params ucp = uc.getParameters();
  for (int i = 0; i < 3; i++) {
    std::cout << ucp.Len(i) << "  " << ucp.Ang(i) << std::endl;
  }
  return 0;
}
