#ifndef SCITBX_ARRAY_FAMILY_SLICE_H
#define SCITBX_ARRAY_FAMILY_SLICE_H

#include <scitbx/array_family/accessors/flex_grid.h>

namespace scitbx { namespace af {

  struct slice {
    std::size_t start, stop, step;

    slice(std::size_t start, std::size_t stop, std::size_t step)
    : start(start), stop(stop), step(step)
    {}

    slice(std::size_t start, std::size_t stop)
    : start(start), stop(stop), step(1)
    {}

  };

  namespace detail {

    template <typename T>
    typename af::versa<T, af::flex_grid<> >::iterator &
    copy_slice_detail(
      af::const_ref<T, af::flex_grid<> > const &self,
      typename af::const_ref<T, af::flex_grid<> >::const_iterator &s,
      typename af::versa<T, af::flex_grid<> >::iterator &r,
      af::small<slice, 10> slices,
      unsigned curr_dim,
      bool copy)
    {
      const unsigned nd = self.accessor().nd();
      slice sl = slices[curr_dim];
      if (curr_dim+1 == nd) {
        if (copy) {
          r = std::copy(s+sl.start, s+sl.stop, r);
        }
        s += self.accessor().all()[curr_dim];
        return r;
      }
      for (unsigned i=0; i<self.accessor().all()[curr_dim]; i++) {
        bool copy_ = copy;
        if (copy_) {
          if (i < sl.start) copy_ = false;
          else if (i >= sl.stop) copy_ = false;
        }
        r = copy_slice_detail(self, s, r, slices, curr_dim+1, copy_);
      }
      return r;
    }

  } // namespace detail

  template <typename T>
  af::versa<T, af::flex_grid<> >
  copy_slice(
    af::const_ref<T, af::flex_grid<> > const &self,
    af::small<slice, 10> slices)
  {
    SCITBX_ASSERT(self.accessor().nd() == slices.size())
      (self.accessor().nd())(slices.size());
    af::small<long, 10> self_dim = self.accessor().all();
    af::small<long, 10> result_dim;
    for (int i=0;i<self.accessor().nd();i++) {
      result_dim.push_back(slices[i].stop-slices[i].start);
    }
    af::versa<T, af::flex_grid<> > result( (af::flex_grid<>(result_dim)) );
    result.resize(af::flex_grid<>(result_dim));
    typename af::versa<T, af::flex_grid<> >::iterator r = result.begin();
    typename af::versa<T, af::flex_grid<> >::const_iterator s = self.begin();
    detail::copy_slice_detail(self, s, r, slices, 0, true);
    return result;
  }

}}

#endif // SCITBX_ARRAY_FAMILY_SLICE_H
