#ifndef CCTBX_SGTBX_RT_MX_HASH_H
#define CCTBX_SGTBX_RT_MX_HASH_H

#include <cctbx/sgtbx/rot_mx_hash.h>
#include <cctbx/sgtbx/tr_vec_hash.h>

namespace cctbx { namespace sgtbx {

  inline
  std::size_t hash_value(rt_mx const &rt) {
    std::size_t result = 0;
    boost::hash_combine(result, rt.t());
    boost::hash_combine(result, rt.r());
    return result;
  }

}}

#endif // GUARD
