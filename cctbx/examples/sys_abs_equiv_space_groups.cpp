#include <cctbx/sgtbx/sys_abs_equiv.h>
#include <iostream>

int
main()
{
  using std::cout;
  using std::endl;
  unsigned const** tab = cctbx::sgtbx::sys_abs_equiv::space_group_numbers;
  // sys_abs_equiv::space_group_numbers has some NULL pointers (0)
  // newer versions of clang do not like dereferencing NULL pointers
  unsigned const zero = 0;
#ifdef __clang__
  if (tab[0] == 0) {
    cout << zero << endl;
  }
#else
  cout << tab[0] << endl;
#endif
  cout << tab[1][1] << endl;
  cout << tab[2][1] << endl;
  for(unsigned i=1;i<=tab[75][0];i++) {
    cout << " " << tab[75][i];
  }
  cout << endl;
#ifdef __clang__
  if (tab[230] == 0) {
    cout << zero << endl; // XXX see above
  }
#else
  cout << tab[230] << endl;
#endif
  return 0;
}
