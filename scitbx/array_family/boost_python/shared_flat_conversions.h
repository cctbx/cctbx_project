/// Conversions flattening and unflattening e.g. shared< tiny<T,N> >.

#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_SHARED_FLAT_CONVERSIONS_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_SHARED_FLAT_CONVERSIONS_H

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/boost_python/utils.h>
#include <boost/python/object.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/to_python_converter.hpp>

namespace scitbx { namespace af { namespace boost_python {

  /// Conversion of e.g. shared< tiny<T,N> > to shared<T>
  /** ElementType is the like of tiny<T,N> */
  template <class ElementType>
  struct structured_to_flat_shared_conversion
  {
    typedef af::shared<ElementType> source_t;
    typedef af::shared<typename ElementType::value_type> target_t;

    static PyObject *convert(source_t const &a) {
      using namespace boost::python;
      target_t result(const_cast<source_t &>(a).handle());
      return incref(object(result).ptr());
    }

    static PyTypeObject const *get_pytype() {
      using namespace boost::python;
      return converter::registered<source_t>::converters.to_python_target_type();
    }

    /// Register conversions
    structured_to_flat_shared_conversion() {
      using namespace boost::python;
      to_python_converter<source_t,
                          structured_to_flat_shared_conversion
                          #ifdef BOOST_PYTHON_SUPPORTS_PY_SIGNATURES
                          , true
                          #endif
                          >();
    }
  };

  /// Conversion of Python flex.T to e.g. const_ref< tiny<T,N> >
  /** ElementType is the like of tiny<T,N> and RefType shall be
      either const_ref<ElementType> or ref<ElementType>.
   */
  template <class ElementType, class RefType>
  struct flat_shared_to_structured_ref_conversion
  {
    typedef af::versa<typename ElementType::value_type, flex_grid<> > source_t;
    typedef RefType target_t;

    static void *convertible(PyObject *obj) {
      using namespace boost::python;
      handle<> hdl(borrowed(obj));
      object py_obj(hdl);
      extract<source_t const &> proxy(py_obj);
      if (!proxy.check()) return 0;
      source_t const &a = proxy();
      if (a.size() % ElementType::size()) return 0;
      return obj;
    }

    static
    void
    construct(PyObject *obj,
              boost::python::converter::rvalue_from_python_stage1_data *data)
    {
      using namespace boost::python;
      handle<> hdl(borrowed(obj));
      object py_obj(hdl);
      source_t &a = extract<source_t &>(py_obj)();
      void
      *storage = ((converter::rvalue_from_python_storage<target_t> *)data
                  )->storage.bytes;
      new (storage) target_t((ElementType *)(&a[0]),
                             a.size()/ElementType::size());
      data->convertible = storage;
    }

    /// Register conversions
    flat_shared_to_structured_ref_conversion() {
      using namespace boost::python;
      converter::registry::push_back(&convertible, &construct,
                                     type_id<target_t>()
                                     #ifdef BOOST_PYTHON_SUPPORTS_PY_SIGNATURES
                                     , &converter::expected_pytype_for_arg<
                                                        source_t>::get_pytype
                                     #endif
                                     );
    }
  };

  /** Conversions to flatten e.g. C++ shared< tiny<T,N> > to Python flex
      and Python flex to C++ ref< tiny<T,N> > and const_ref< tiny<T,N> >.

   Thanks to these conversions, a shared< tiny<T,N> > can be passed from C++
   to Python and back without the need of a Boost.Python wrapper for that array
   type, thus reducing code bloat while retaining convenience.
   */
  template <class ElementType>
  struct flat_shared_conversions
  {
    flat_shared_conversions() {
      structured_to_flat_shared_conversion<ElementType>();
      flat_shared_to_structured_ref_conversion<ElementType,
                                               ref<ElementType> >();
      flat_shared_to_structured_ref_conversion<ElementType,
                                               const_ref<ElementType> >();
    }
  };

}}}

#endif // GUARD
