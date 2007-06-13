#include <scitbx/error.h>
#import <string>

int main() {
  double x = 1.1;
  int n = 1;
  SCITBX_ASSERT(x*x*x < n)(x)(n);
};

