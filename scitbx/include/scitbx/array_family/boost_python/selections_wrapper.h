#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_SELECTIONS_WRAPPER_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_SELECTIONS_WRAPPER_H

#include <boost/python/args.hpp>
#include <scitbx/array_family/selections.h>

namespace scitbx { namespace af { namespace boost_python {

  template <typename ElementType, typename SelfType>
  struct select_wrappers
  {
    static shared<ElementType>
    with_flags(
      SelfType const& self,
      const_ref<bool> const& flags)
    {
      return select(self.const_ref().as_1d(), flags);
    }

    static shared<ElementType>
    with_indices(
      SelfType const& self,
      const_ref<std::size_t> const& indices)
    {
      return select(self.const_ref().as_1d(), indices, false);
    }

    static shared<ElementType>
    with_indices_reverse(
      SelfType const& self,
      const_ref<std::size_t> const& indices,
      bool reverse)
    {
      return select(self.const_ref().as_1d(), indices, reverse);
    }

    template <typename ArrayWrapper>
    static void
    wrap(ArrayWrapper& aw)
    {
      using namespace boost::python;
      aw.def("select", with_flags, (arg_("self"), arg_("flags")))
        .def("select", with_indices, (arg_("self"), arg_("indices")))
        .def("select", with_indices_reverse, (
          arg_("self"), arg_("indices"), arg_("reverse")))
      ;
    }
  };

}}} // namespace scitbx::af::boost_python

#endif // SCITBX_ARRAY_FAMILY_BOOST_PYTHON_SELECTIONS_WRAPPER_H
