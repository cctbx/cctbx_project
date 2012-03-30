#ifndef SCITBX_ARRAY_FAMILY_SELECTIONS_H
#define SCITBX_ARRAY_FAMILY_SELECTIONS_H

#include <scitbx/array_family/shared.h>
#include <scitbx/error.h>

namespace scitbx { namespace af {

  template <
    typename ElementType,
    typename UnsignedType>
  shared<ElementType>
  select(
    const_ref<ElementType> const& self,
    const_ref<UnsignedType> const& indices,
    bool reverse=false)
  {
    if (!reverse) {
      shared<ElementType> result((reserve(indices.size())));
      for(std::size_t i=0;i<indices.size();i++) {
        SCITBX_ASSERT(indices[i] < self.size());
        result.push_back(self[indices[i]]);
      }
      return result;
    }
    SCITBX_ASSERT(indices.size() == self.size());
    shared<ElementType> result;
    if (self.size()) {
      result.resize(self.size(), self[0]); // avoid requirement that element
      for(std::size_t i=1;i<self.size();i++) {   // is default constructible
        SCITBX_ASSERT(indices[i] < self.size());
        result[indices[i]] = self[i];
      }
    }
    return result;
  }

  template <typename ElementType>
  shared<ElementType>
  select(
    const_ref<ElementType> const& self,
    const_ref<bool> const& flags)
  {
    SCITBX_ASSERT(flags.size() == self.size());
    std::size_t n = 0;
    for(std::size_t i=0;i<flags.size();i++) if (flags[i]) n++;
    shared<ElementType> result((reserve(n)));
    for(std::size_t i=0;i<flags.size();i++) {
      if (flags[i]) result.push_back(self[i]);
    }
    return result;
  }

  template <typename IntType, typename IntTypeOther>
  struct intersection_with_tracking
  {
    shared<IntType> matching_elements;
    shared<std::size_t> self_i_seqs;
    shared<std::size_t> other_i_seqs;

    intersection_with_tracking() {}

    intersection_with_tracking(
      const_ref<IntType> const& self,
      const_ref<IntTypeOther> const& other,
      bool track_matching_elements,
      bool track_i_seqs)
    {
      static const char* err1_dup = "the first array has duplicate elements.";
      static const char* err1_not = "the first array is not sorted.";
      static const char* err2_dup = "the second array has duplicate elements.";
      static const char* err2_not = "the second array is not sorted.";
      if (self.size() > 0 && other.size() > 0) {
        IntType si = self[0];
        std::size_t i = 1;
        IntTypeOther oj = other[0];
        std::size_t j = 1;
        while (true) {
          while (si < oj) {
            if (i == self.size()) return;
            if (si == self[i]) throw SCITBX_ERROR(err1_dup);
            if (si > self[i]) throw SCITBX_ERROR(err1_not);
            si = self[i];
            i++;
          }
          while (oj < si) {
            if (j == other.size()) return;
            if (oj == other[j]) throw SCITBX_ERROR(err2_dup);
            if (oj > other[j]) throw SCITBX_ERROR(err2_not);
            oj = other[j];
            j++;
          }
          if (oj == si) {
            if (track_matching_elements) {
              matching_elements.push_back(si);
            }
            if (track_i_seqs) {
              self_i_seqs.push_back(i-1);
              other_i_seqs.push_back(j-1);
            }
            if (i == self.size()) break;
            if (si == self[i]) throw SCITBX_ERROR(err1_dup);
            if (si > self[i]) throw SCITBX_ERROR(err1_not);
            si = self[i];
            i++;
          }
        }
      }
    }
  };

  template <typename IntType, typename IntTypeOther>
  shared<IntType>
  intersection(
    const_ref<IntType> const& self,
    const_ref<IntTypeOther> const& other)
  {
    return intersection_with_tracking<IntType, IntTypeOther>(
      self,
      other,
      /*track_matching_elements*/ true,
      /*track_i_seqs*/ false).matching_elements;
  }

  template <typename IntType>
  struct reindexing_given_bool_selection
  {
    IntType n_selected;
    shared<IntType> array;

    reindexing_given_bool_selection() {}

    reindexing_given_bool_selection(
      const_ref<bool> const& selection)
    :
      n_selected(0),
      array(selection.size(), selection.size())
    {
      IntType* array_begin = array.begin();
      for(std::size_t i=0;i<selection.size();i++) {
        if (!selection[i]) continue;
        array_begin[i] = n_selected++;
      }
    }
  };

  template <typename IntType>
  shared<IntType>
  reindexing_array(
    std::size_t selectee_size,
    const_ref<IntType> const& iselection)
  {
    shared<IntType> result(selectee_size, selectee_size);
    IntType* result_begin = result.begin();
    for(std::size_t i=0;i<iselection.size();i++) {
      SCITBX_ASSERT(iselection[i] < selectee_size);
      result_begin[iselection[i]] = i;
    }
    return result;
  }

  template <typename MapType>
  shared<MapType>
  array_of_map_proxy_select(
    const_ref<MapType> const& self,
    const_ref<std::size_t> const& iselection)
  {
    std::size_t selectee_size = self.size();
    shared<std::size_t> reindexing_array = af::reindexing_array(
      selectee_size, iselection);
    std::size_t* reindexing_array_begin = reindexing_array.begin();
    shared<MapType> result((reserve(iselection.size())));
    for(std::size_t i=0;i<iselection.size();i++) {
      result.push_back(MapType());
      MapType& new_map = result.back();
      MapType const& old_map = self[iselection[i]];
      for(typename MapType::const_iterator
            old_map_i=old_map.begin();
            old_map_i!=old_map.end();
            old_map_i++) {
        SCITBX_ASSERT(old_map_i->first < selectee_size);
        std::size_t j = reindexing_array_begin[old_map_i->first];
        if (j != selectee_size) {
          new_map[static_cast<typename MapType::key_type>(j)]
            = old_map_i->second;
        }
      }
    }
    return result;
  }

  template <typename MapType>
  shared<MapType>
  array_of_map_proxy_remove(
    const_ref<MapType> const& self,
    const_ref<bool> const& selection)
  {
    SCITBX_ASSERT(selection.size() == self.size());
    shared<MapType> result;
    for(std::size_t i=0;i<self.size();i++) {
      if (!selection[i]) {
        result.push_back(self[i]);
      }
      else {
        result.push_back(MapType());
        MapType& new_map = result.back();
        MapType const& old_map = self[i];
        for(typename MapType::const_iterator
              old_map_i=old_map.begin();
              old_map_i!=old_map.end();
              old_map_i++) {
          SCITBX_ASSERT(old_map_i->first < self.size());
          if (!selection[old_map_i->first]) {
            new_map[old_map_i->first] = old_map_i->second;
          }
        }
      }
    }
    return result;
  }

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_SELECTIONS_H
