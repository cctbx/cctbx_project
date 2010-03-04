#ifndef CCTBX_SGTBX_TR_GROUP_HASH_H
#define CCTBX_SGTBX_TR_GROUP_HASH_H

#include <cctbx/sgtbx/tr_group.h>

namespace cctbx { namespace sgtbx {

  inline
  std::size_t hash_value(tr_group const &trg) {
    std::size_t result = 0;
    boost::hash_combine(result, trg.t_den());
    boost::hash_range(result, trg.elems().begin(), trg.elems().end());
    return result;
  }

}}

#endif // GUARD
