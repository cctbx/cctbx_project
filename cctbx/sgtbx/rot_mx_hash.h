#ifndef CCTBX_SGTBX_ROT_MX_HASH_H
#define CCTBX_SGTBX_ROT_MX_HASH_H

#include <boost/functional/hash.hpp>

namespace cctbx { namespace sgtbx {

  inline
  std::size_t hash_value(rot_mx const &r) {
    std::size_t result = r.den();
    boost::hash_range(result, r.num().begin(), r.num().end());
    return result;
  }

}}

#endif // GUARD
