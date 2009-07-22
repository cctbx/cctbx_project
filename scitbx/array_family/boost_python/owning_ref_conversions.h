/// Convertion of the owning references defined in
/// scitbx/array_family/owning_ref.h to Python

#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_OWNING_REF_CONVERSIONS_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_OWNING_REF_CONVERSIONS_H

#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/owning_ref.h>

#include <boost/python/object.hpp>
#include <boost/python/extract.hpp>

namespace scitbx { namespace af { namespace boost_python {

  template <class T>
  struct ref_owning_shared_to_flex
  {
    static PyObject *convert(ref_owning_shared<T> const &a) {
      using namespace boost::python;
      versa<T, flex_grid<> > result(a.array(), flex_grid<>(a.size()));
      return incref(object(result).ptr());
    }

    static const PyTypeObject *get_pytype() {
      using namespace boost::python;
      return converter::registered< versa<T, flex_grid<> > >
              ::converters.to_python_target_type();
    }
  };

  template <class T, class A>
  struct ref_owning_versa_to_flex
  {
    static PyObject *convert(ref_owning_versa<T, A> const &a) {
      using namespace boost::python;
      versa<T, flex_grid<> > result(a.array(), a.accessor().as_flex_grid());
      return incref(object(a.array()).ptr());
    }

    static const PyTypeObject *get_pytype() {
      using namespace boost::python;
      return converter::registered< versa<T, flex_grid<> > >
              ::converters.to_python_target_type();
    }
  };

  template <class T>
  struct ref_owning_shared_conversions
  {
    ref_owning_shared_conversions() {
      using namespace boost::python;
      to_python_converter<ref_owning_shared<T>,
                          ref_owning_shared_to_flex<T>
                          #ifdef BOOST_PYTHON_SUPPORTS_PY_SIGNATURES
                            , true
                          #endif
                          >();
    }
  };

  template <class T, class A>
  struct ref_owning_versa_conversions
  {
    ref_owning_versa_conversions() {
      using namespace boost::python;
      to_python_converter<ref_owning_versa<T, A>,
                          ref_owning_versa_to_flex<T, A>
                          #ifdef BOOST_PYTHON_SUPPORTS_PY_SIGNATURES
                            , true
                          #endif
                          >();
    }
  };

}}}

#endif // GUARD
