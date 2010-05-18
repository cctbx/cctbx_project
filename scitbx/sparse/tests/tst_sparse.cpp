#include <scitbx/sparse/vector.h>
#include <iostream>

namespace scitbx { namespace sparse {

  /// Test the efficient assignment and augmented assignment available in C++
  void exercise_vector_element_assignment() {
    {
      vector<double> v;
      SCITBX_ASSERT(!v.is_compact());
      v.compact();
      SCITBX_ASSERT(v.is_compact());
      SCITBX_ASSERT(v.size() == 0);
    }
    {
      vector<double> v(4);
      SCITBX_ASSERT(!v.is_compact());
      for (int i=0; i<4; ++i) SCITBX_ASSERT(v[i] == 0)(v[i]);
      v.compact();
      SCITBX_ASSERT(v.is_compact());
      SCITBX_ASSERT(v.size() == 4)(v.size());
    }
    {
      unsigned n = 5;
      vector<double> v(n);
      v[1] += 1e-20;
      v[2] += 2.;
      v.compact();
      af::shared<double> w(n);
      w[2] = 2; w[1] = 1e-20;
      for (vector<double>::iterator p=v.begin(); p != v.end(); ++p) {
        SCITBX_ASSERT(*p == w[p.index()])(p.index())(*p);
      }
    }
    {
      vector<double> v(6);
      v[4] += 1.;
      v[5] += 2.;
      v[4] += 3.;
      v[5] = 4.;
      v[3] = 5.;
      v[5] += 6.;
      SCITBX_ASSERT(!v.is_compact());
      v.compact();
      SCITBX_ASSERT(v.is_compact());
      SCITBX_ASSERT(v.size() == 6);
      af::shared<double> w(6);
      w[4] = 4.; w[5] = 10.; w[3] = 5.;
      for (vector<double>::iterator p=v.begin(); p != v.end(); ++p) {
        SCITBX_ASSERT(*p == w[p.index()])(p.index())(*p);
      }
      v[1] = 0;
      SCITBX_ASSERT(!v.is_compact());
    }
    {
      vector<double> v;
      v[2] += 1.;
      v[3] += 2.;
      v[5] += 3.;
      v[3] += 4.;
      v[2] += 5.;
      v[5] += 6.;
      SCITBX_ASSERT(!v.is_compact());
      v.compact();
      SCITBX_ASSERT(v.is_compact());
      SCITBX_ASSERT(v.size() == 6);
      af::shared<double> w(6);
      w[2] = 6.; w[3] = 6.; w[5] = 9.;
      for (vector<double>::iterator p=v.begin(); p != v.end(); ++p) {
        SCITBX_ASSERT(*p == w[p.index()])(p.index())(*p);
      }
    }
    {
      vector<double> v;
      v[1] = 1.;
      v[4] = 2.;
      v[6] = 3.;
      SCITBX_ASSERT(!v.is_compact());
      v.compact();
      SCITBX_ASSERT(v.is_compact());
      SCITBX_ASSERT(v.size() == 7);
      af::shared<double> w(7);
      w[1] = 1.; w[4] = 2.; w[6] = 3.;
      for (vector<double>::iterator p=v.begin(); p != v.end(); ++p) {
        SCITBX_ASSERT(*p == w[p.index()])(p.index())(*p);
      }
    }
  }
}}

int main() {
  scitbx::sparse::exercise_vector_element_assignment();
  std::cout << "OK\n" << std::endl;
  return 0;
}
