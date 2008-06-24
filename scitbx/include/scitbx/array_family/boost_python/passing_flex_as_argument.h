#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_PASSING_FLEX_AS_ARGUMENT_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_PASSING_FLEX_AS_ARGUMENT_H

#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/versa.h>

namespace scitbx { namespace af {

template<class ElementType>
class flex_1d_argument
{
  public:
    typedef af::versa<ElementType, af::flex_grid<> > array_type;
    typedef typename array_type::base_array_type base_array_type;

    flex_1d_argument(array_type &array)
      : a(array), b(array.as_base_array())
    {
      SCITBX_ASSERT(array.accessor().nd() == 1 && array.accessor().is_0_based())
                   (array.accessor().nd());
    }

    ~flex_1d_argument() { a.resize(af::flex_grid<>(b.size())); }

    base_array_type *operator->() { return &b; }

    ElementType &operator[](std::size_t i) { return b[i]; }

  private:
    array_type &a;
    base_array_type b;
};

}} // namespace scitbx::af

#endif // GUARD
