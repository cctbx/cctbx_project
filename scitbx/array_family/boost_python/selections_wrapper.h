#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_SELECTIONS_WRAPPER_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_SELECTIONS_WRAPPER_H

#include <boost/python/args.hpp>
#include <boost_adaptbx/boost_python_type_id_eq.h>
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

    template <typename UnsignedType>
    static shared<ElementType>
    with_indices(
      SelfType const& self,
      const_ref<UnsignedType> const& indices,
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
        .def("select", with_indices<unsigned>, (
          arg_("self"), arg_("indices"), arg_("reverse")=false))
#if !defined(BOOST_PYTHON_TYPE_ID_UNSIGNED_EQ_SIZE_T)
        .def("select", with_indices<std::size_t>, (
          arg_("self"), arg_("indices"), arg_("reverse")=false))
#endif
      ;
    }
  };

}}} // namespace scitbx::af::boost_python

#endif // SCITBX_ARRAY_FAMILY_BOOST_PYTHON_SELECTIONS_WRAPPER_H
