#ifndef CCTBX_SGTBX_SPACE_GROUP_HASH_H
#define CCTBX_SGTBX_SPACE_GROUP_HASH_H

#include <cctbx/sgtbx/space_group.h>
#include <cctbx/sgtbx/rt_mx_hash.h>
#include <cctbx/sgtbx/tr_group_hash.h>

namespace cctbx { namespace sgtbx {

  inline
  std::size_t hash_value(space_group const &sg) {
    if (!sg.is_tidy()) {
      throw std::runtime_error("Make space group tidy before hashing it");
    }
    std::size_t result = 0;
    boost::hash_combine(result, sg.r_den());
    boost::hash_combine(result, sg.t_den());
    boost::hash_combine(result, sg.is_centric());
    if (sg.is_centric()) boost::hash_combine(result, sg.inv_t());
    boost::hash_combine(result, sg.ltr());
    boost::hash_range(result, sg.smx().begin(), sg.smx().end());
    return result;
  }

}}

#endif // GUARD
