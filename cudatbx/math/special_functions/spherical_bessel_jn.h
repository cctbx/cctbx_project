#ifndef CUDATBX_SPHERICAL_BESSEL_JN_H
#define CUDATBX_SPHERICAL_BESSEL_JN_H

#include <scitbx/array_family/shared.h>

namespace cudatbx {
namespace math {
namespace special_functions {

  scitbx::af::shared<double> cuda_spherical_bessel_jn
    (const int&, const scitbx::af::const_ref<double>&, const int&);

}
}
}
#endif // CUDATBX_SPHERICAL_BESSEL_JN_H
