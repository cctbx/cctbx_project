#include <scitbx/boost_python/slice.h>
#include <scitbx/error.h>

namespace scitbx { namespace boost_python {

  adapted_slice::adapted_slice(slice const& sl, std::size_t sz)
  :
    step(1),
    size(0)
  {
    long signed_sz = static_cast<long>(sz);
    if (sl.step.is_valid) {
      step = sl.step.value;
    }
    if (!sl.start.is_valid) {
      start = step < 0 ? signed_sz-1 : 0;
    }
    else {
      start = sl.start.value;
      if (start < 0) start += signed_sz;
    }
    if (!sl.stop.is_valid) {
      stop = step < 0 ? -1 : signed_sz;
    }
    else {
      stop = sl.stop.value;
      if (stop < 0) stop += signed_sz;
    }
    if (start > signed_sz-1) start = signed_sz;
    if (start < 0) start = 0;
    if      (stop < -1) stop = -1;
    else if (stop > signed_sz) stop = signed_sz;
    SCITBX_ASSERT(step != 0 || stop == start);
    long signed_size = stop - start + step;
    if      (step < 0) signed_size++;
    else if (step > 0) signed_size--;
    else {
      signed_size = 0;
      step = 1;
    }
    signed_size /= step;
    if (signed_size < 0) signed_size = 0;
    size = static_cast<std::size_t>(signed_size);
    stop = start + step * size;
  }

}} // namespace scitbx::boost_python
