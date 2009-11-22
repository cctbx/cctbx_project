#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_SELECTIONS_WRAPPER_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_SELECTIONS_WRAPPER_H

#include <boost/python/args.hpp>
#include <boost_adaptbx/type_id_eq.h>
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

    // gcc 3.2 fails with template <typename UnsignedType>
    static shared<ElementType>
    with_indices_unsigned(
      SelfType const& self,
      const_ref<unsigned> const& indices,
      bool reverse)
    {
      return select(self.const_ref().as_1d(), indices, reverse);
    }

    static shared<ElementType>
    with_indices_size_t(
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
      using boost::python::arg;
      aw.def("select", with_flags, (arg("self"), arg("flags")))
        .def("select", with_indices_unsigned, (
          arg("self"), arg("indices"), arg("reverse")=false))
#if !defined(BOOST_ADAPTBX_TYPE_ID_SIZE_T_EQ_UNSIGNED)
        .def("select", with_indices_size_t, (
          arg("self"), arg("indices"), arg("reverse")=false))
#endif
      ;
    }
  };

}}} // namespace scitbx::af::boost_python

#endif // SCITBX_ARRAY_FAMILY_BOOST_PYTHON_SELECTIONS_WRAPPER_H
