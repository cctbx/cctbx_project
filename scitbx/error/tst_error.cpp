#include <scitbx/error.h>
#include <boost_adaptbx/error_utils.h>
#include <iostream>

int
main()
{
  {
    double x = 1.1;
    int n = 1;
    int j = 2;
    bool have_message = false;
    try {
      SCITBX_ASSERT(x*x*x < n)(x)(n)(j);
    }
    catch (scitbx::error const& e) {
      std::cout << e.what() << std::endl;
      have_message = true;
    }
    SCITBX_ASSERT(have_message);
  }
  {
    bool have_message = false;
    try {
      throw BOOST_ADAPTBX_UNREACHABLE_ERROR();
    }
    catch (std::runtime_error const& e) {
      std::cout << e.what() << std::endl;
      have_message = true;
    }
    SCITBX_ASSERT(have_message);
  }
  return 0;
}
