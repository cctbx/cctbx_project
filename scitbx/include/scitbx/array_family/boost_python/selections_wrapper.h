#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_SELECTIONS_WRAPPER_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_SELECTIONS_WRAPPER_H

#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>
#include <scitbx/array_family/selections.h>

namespace scitbx { namespace af { namespace boost_python {

  template <typename ElementType>
  struct select_wrappers
  {
    BOOST_PYTHON_FUNCTION_OVERLOADS(select_overloads, select, 2, 3)

    template <typename ArrayWrapper>
    static void
    wrap(ArrayWrapper& aw)
    {
      using namespace boost::python;
      aw.def("select",
          (shared<ElementType>(*)(
            const_ref<ElementType> const&,
            const_ref<bool> const&)) select, (
              arg_("self"), arg_("flags")))
        .def("select",
          (shared<ElementType>(*)(
            const_ref<ElementType> const&,
            const_ref<std::size_t> const&,
            bool)) select, select_overloads((
              arg_("self"), arg_("indices"), arg_("reverse")=false)))
      ;
    }
  };

}}} // namespace scitbx::af::boost_python

#endif // SCITBX_ARRAY_FAMILY_BOOST_PYTHON_SELECTIONS_WRAPPER_H
