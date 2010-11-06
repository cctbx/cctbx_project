#ifndef CCTBX_SGTBX_TR_VEC_HASH_H
#define CCTBX_SGTBX_TR_VEC_HASH_H

#include <cctbx/sgtbx/tr_vec.h>
#include <boost/functional/hash.hpp>

namespace cctbx { namespace sgtbx {

  inline
  std::size_t hash_value(tr_vec const &v) {
    std::size_t result = v.den();
    boost::hash_range(result, v.num().begin(), v.num().end());
    return result;
  }

}}

#endif // GUARD
