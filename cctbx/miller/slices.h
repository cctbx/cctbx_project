#ifndef CCTBX_MILLER_SLICES_H
#define CCTBX_MILLER_SLICES_H

#include <cctbx/miller.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/ref.h>

namespace cctbx { namespace miller {

  af::shared<bool>
  simple_slice (
    af::const_ref<index<int> > const& indices,
    unsigned slice_axis,
    int slice_index)
  {
    CCTBX_ASSERT((slice_axis >= 0) && (slice_axis < 3));
    af::shared<bool> selection(indices.size(), false);
    for (std::size_t i_seq = 0; i_seq < indices.size(); i_seq++) {
      if (indices[i_seq][slice_axis] == slice_index) {
        selection[i_seq] = true;
      }
    }
    return selection;
  }

  af::shared<bool>
  multi_slice (
    af::const_ref<index<int> > const& indices,
    unsigned slice_axis,
    int slice_start,
    int slice_end)
  {
    CCTBX_ASSERT((slice_axis >= 0) && (slice_axis < 3));
    CCTBX_ASSERT((slice_start <= slice_end));
    af::shared<bool> selection(indices.size(), false);
    for (std::size_t i_seq = 0; i_seq < indices.size(); i_seq++) {
      if ((indices[i_seq][slice_axis] >= slice_start) &&
          (indices[i_seq][slice_axis] <= slice_end)) {
        selection[i_seq] = true;
      }
    }
    return selection;
  }

}}
#endif
