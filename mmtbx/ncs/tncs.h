#ifndef MMTBX_TNCS_H
#define MMTBX_TNCS_H

#include <cctbx/import_scitbx_af.h>
#include <boost/python/list.hpp>
#include <scitbx/sym_mat3.h>

using namespace std;
namespace mmtbx { namespace ncs {
namespace af=scitbx::af;
using scitbx::vec3;
using scitbx::mat3;

template <typename FloatType=double>
class pair
{
  public:
    scitbx::mat3<FloatType> r;
    scitbx::vec3<FloatType> t;
    FloatType radius;
    FloatType weight;
    pair() {
     r.fill(0.0);
     t.fill(0.0);
     radius=0;
     weight=0;
    }
    pair(scitbx::mat3<FloatType> const& r_,
         scitbx::vec3<FloatType> const& t_,
         FloatType radius_,
         FloatType weight_)
    :
      r(r_),t(t_),radius(radius_),weight(weight_)
    {}
};

}} // namespace mmtbx::ncs

#endif // MMTBX_TNCS_H
