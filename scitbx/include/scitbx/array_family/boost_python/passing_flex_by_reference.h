#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_PASSING_FLEX_AS_ARGUMENT_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_PASSING_FLEX_AS_ARGUMENT_H

#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/versa.h>

namespace scitbx { namespace af {

template<class ElementType>
class flex_1d_argument
  : public af::versa<ElementType, af::flex_grid<> >::base_array_type
{
  public:
    typedef af::versa<ElementType, af::flex_grid<> > flex_array_type;
    typedef typename flex_array_type::base_array_type base_type;

    flex_1d_argument(flex_array_type &array)
      : a(array), base_type(array.as_base_array())
    {
      SCITBX_ASSERT(array.accessor().nd() == 1 && array.accessor().is_0_based())
                   (array.accessor().nd());
    }

    ~flex_1d_argument() { a.resize(af::flex_grid<>(this->size())); }

  private:
    flex_array_type &a;
};

}} // namespace scitbx::af

#endif // GUARD
