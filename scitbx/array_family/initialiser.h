#ifndef SCITBX_ARRAY_FAMILY_INITIALISER_H
#define SCITBX_ARRAY_FAMILY_INITIALISER_H

#include <scitbx/error.h>

namespace scitbx { namespace af {

/// Syntactic sugar to initialise arrays
/** Synopsis:

    af::versa<double, af::c_grid<2> > m(af::c_grid<2>(3,2));
    af::init(m) = 1, 2,
                  3, 4,
                  5, 6;

    Rationale: especially useful to write test cases.
*/
template <class IndexedType>
struct initialiser
{
  typedef typename IndexedType::value_type value_type;

  IndexedType &a;
  std::size_t i;

  initialiser(IndexedType &indexed) : a(indexed), i(0) {}

  initialiser &operator=(value_type const &x) {
    SCITBX_ASSERT(i < a.size());
    a[i++] = x;
    return *this;
  }

  initialiser &operator,(value_type const &x) {
    SCITBX_ASSERT(i < a.size());
    a[i++] = x;
    return *this;
  }
};


template <class IndexedType>
initialiser<IndexedType> init(IndexedType &indexed) {
  return initialiser<IndexedType>(indexed);
}

}} // scitbx::af

#endif // GUARD
