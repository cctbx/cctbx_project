#include <iostream>
#include <iotbx/shelx/hklf.h>
#include <scitbx/sym_mat3.h>
#include <scitbx/mat3.h>
#include <scitbx/vec3.h>
#include <scitbx/math/approx_equal.h>
#include <iotbx/error.h>

class test_case1 {
  public:
    typedef scitbx::vec3<double> vec3_t;
    typedef scitbx::mat3<double> mat3_t;
    typedef scitbx::sym_mat3<double> sym_mat3_t;

    void run() {
      std::cout<<"########################################"<<std::endl;
      vec3_t a(1,2,3);
      vec3_t b(4,5,6);
      sym_mat3_t g(1,1,1,0,0,0);
      vec3_t v(-3,6,-3);

      vec3_t out = iotbx::shelx::hklf_reader::make_perpendicular(a, b, g);
      v = v.normalize();
      for (int i=0; i<3; i++) {
        IOTBX_ASSERT(scitbx::math::approx_equal_absolutely<double>(1e-12)(v[i], out[i]));
      }
      IOTBX_ASSERT(scitbx::math::approx_equal_absolutely<double>(1e-12)(out.length(), 1));


      g = sym_mat3_t(10,11,12,1,2,3);
      out = iotbx::shelx::hklf_reader::make_perpendicular(a, b, g);
      IOTBX_ASSERT(scitbx::math::approx_equal_absolutely<double>(1e-12)(0, a*g*out));
      IOTBX_ASSERT(scitbx::math::approx_equal_absolutely<double>(1e-12)(0, b*g*out));

    }
};



int main()
{
  using namespace iotbx::shelx;


  test_case1 t;
  t.run();

  return 0;
}
