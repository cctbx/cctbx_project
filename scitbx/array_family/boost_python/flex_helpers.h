#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_FLEX_HELPERS_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_FLEX_HELPERS_H

#include <boost/python/args.hpp>

namespace scitbx { namespace af { namespace boost_python {

  template <typename ElementType, typename UnsignedType>
  boost::python::object
  add_selected_unsigned_a(
    boost::python::object const& self,
    af::const_ref<UnsignedType> const& indices,
    af::const_ref<ElementType> const& values)
  {
    boost::python::extract<af::ref<ElementType> > a_proxy(self);
    af::ref<ElementType> a = a_proxy();
    SCITBX_ASSERT(indices.size() == values.size());
    for(std::size_t i=0;i<indices.size();i++) {
      SCITBX_ASSERT(indices[i] < a.size());
      a[indices[i]] += values[i];
    }
    return self;
  }

  struct approx_equal_helper
  {
    template <typename ElementType, class ConstRefOrElementType>
    static
    bool all_approx_equal(
      const_ref<ElementType> const& self,
      ConstRefOrElementType const& other,
      typename const_ref<ElementType>::amplitude_type tolerance=1.e-6)
    {
      return self.all_approx_equal(other, tolerance);
    }

    template <typename ElementType, class ConstRefOrElementType>
    static
    bool all_approx_equal_relatively(
      const_ref<ElementType> const& self,
      ConstRefOrElementType const& other,
      typename const_ref<ElementType>::amplitude_type tolerance=1.e-6)
    {
      return self.all_approx_equal_relatively(other, tolerance);
    }

    template <class BoostPythonClass>
    static
    void decorate(BoostPythonClass &klass) {
      typedef typename BoostPythonClass::wrapped_type f_t;
      typedef typename f_t::value_type                e_t;

      using boost::python::arg;

      klass
        .def("all_approx_equal",
             all_approx_equal<e_t, const_ref<e_t> >,
             (arg("self"),
              arg("other"),
              arg("tolerance")=1.e-6))
        .def("all_approx_equal",
             all_approx_equal<e_t, e_t>,
             (arg("self"),
              arg("other"),
              arg("tolerance")=1.e-6))
        .def("all_approx_equal_relatively",
             all_approx_equal_relatively<e_t, const_ref<e_t> >,
             (arg("self"),
              arg("other"),
              arg("relative_error")))
        .def("all_approx_equal_relatively",
             all_approx_equal_relatively<e_t, e_t>,
             (arg("self"),
              arg("other"),
              arg("relative_error")))
        ;
    }
  };



}}} // namespace scitbx::af::boost_python



#endif
