#include <cctbx/sgtbx/sys_abs_equiv.h>
#include <iostream>

int
main()
{
  using std::cout;
  using std::endl;
  unsigned const** tab = cctbx::sgtbx::sys_abs_equiv::space_group_numbers;
  cout << (*tab[0]) << endl; // XXX need to dereference when using clang
  cout << tab[1][1] << endl;
  cout << tab[2][1] << endl;
  for(unsigned i=1;i<=tab[75][0];i++) {
    cout << " " << tab[75][i];
  }
  cout << endl;
  cout << (*tab[230]) << endl; // XXX see above
  return 0;
}
