// $Id$

#include <cctbx/uctbx.h>
#include <cctbx/miller.h>

int main(void)
{
  uctbx::UnitCell uc(uctbx::uc_params(2, 3, 4, 81, 82, 83));
  std::cout << uc << std::endl;
  uctbx::uc_params ucp = uc.getParameters();
  for (int i = 0; i < 3; i++) {
    std::cout << ucp.Len(i) << "  " << ucp.Ang(i) << std::endl;
  }
  Miller::Index m1(2,3,4);
  Miller::Index m2(3,4,2);
  std::cout << (m1 < m1) << (m1 > m1) << std::endl;
  std::cout << (m1 < m1) << (m1 > m1) << std::endl;
  std::cout << (m1 == m2) << (m1 != m2) << std::endl;
  std::cout << (m1 == m1) << (m1 != m1) << std::endl;
  return 0;
}
