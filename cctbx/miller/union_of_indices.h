#ifndef CCTBX_MILLER_UNION_OF_INDICES_H
#define CCTBX_MILLER_UNION_OF_INDICES_H

#include <cctbx/miller.h>
#include <set>

namespace cctbx { namespace miller {

  struct union_of_indices_registry
  {
    std::set<index<> > storage;

    union_of_indices_registry() {}

    void
    update(
      af::const_ref<index<> > const& indices)
    {
      for(std::size_t i=0;i<indices.size();i++) {
        storage.insert(indices[i]);
      }
    }

    af::shared<index<> >
    as_array() const
    {
      af::shared<index<> > result;
      result.reserve(storage.size());
      typedef std::set<index<> >::const_iterator it;
      it s_end = storage.end();
      for(it i=storage.begin();i!=s_end;i++) {
        result.push_back(*i);
      }
      return result;
    }
  };

}} // namespace cctbx::miller

#endif // GUARD
