// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Jan 14: Created, based on std_vector_bpl.h (rwgk)
 */

#ifndef CCTBX_SHARED_STORAGE_BPL_H
#define CCTBX_SHARED_STORAGE_BPL_H

#include <cctbx/ndim_bpl.h>
#include <cctbx/shared_storage.h>

namespace cctbx {

  void wrap_shared_storage_handle(boost::python::module_builder& this_module)
  {
    boost::python::class_builder<cctbx::shared_storage_handle_type>
    py_shared_storage_handle(this_module, "shared_storage_handle");

    py_shared_storage_handle.def(boost::python::constructor<>());
  }

  // A wrapper is used to define additional constructors.
  template <typename ValueType>
  struct shared_storage_wrapper : cctbx::shared_storage<ValueType>
  {
    // Tell the compiler how to convert a base class object to
    // this wrapper object.
    shared_storage_wrapper(PyObject*,
                           const cctbx::shared_storage<ValueType>& v)
      : cctbx::shared_storage<ValueType>(v)
    {}
    shared_storage_wrapper(PyObject*)
      : cctbx::shared_storage<ValueType>()
    {}
    shared_storage_wrapper(PyObject* self,
                           std::size_t n)
      : cctbx::shared_storage<ValueType>(n)
    {}
    shared_storage_wrapper(PyObject* self,
                           const cctbx::shared_storage_handle_type& handle,
                           std::size_t n)
      : cctbx::shared_storage<ValueType>(handle, n)
    {}
    shared_storage_wrapper(PyObject* self,
                           boost::python::tuple tuple)
      : cctbx::shared_storage<ValueType>(tuple.size())
    {
      cctbx::shared_storage<ValueType>::iterator v = this->begin();
      for (std::size_t i = 0; i < tuple.size(); i++)
        v[i] = BOOST_PYTHON_CONVERSION::from_python(tuple[i].get(),
                                           boost::python::type<ValueType>());
    }

  };

  void raise_IndexError() {
    PyErr_SetString(PyExc_IndexError, "index out of range");
    throw boost::python::error_already_set();
  }

  template <typename ValueType>
  struct shared_storage_access
  {
    static
    ValueType
    getitem(const cctbx::shared_storage<ValueType>& v,
            std::size_t key)
    {
      if (key >= v.size()) raise_IndexError();
      return v[key];
    }

    static
    void
    setitem(cctbx::shared_storage<ValueType>& v,
            std::size_t key,
            const ValueType& value)
    {
      if (key >= v.size()) raise_IndexError();
      v[key] = value;
    }

    static
    cctbx::shared_storage_handle_type&
    handle(cctbx::shared_storage<ValueType>& v)
    {
      return v.handle();
    }

    static
    void fill(cctbx::shared_storage<ValueType>& v, const ValueType& value)
    {
      std::fill(v.begin(), v.end(), value);
    }

    // Convert to a regular Python tuple.
    static
    boost::python::tuple
    as_tuple(const cctbx::shared_storage<ValueType>& v)
    {
      boost::python::tuple t(v.size());
      for (std::size_t i = 0; i < v.size(); i++) {
        t.set_item(i,
          boost::python::ref(BOOST_PYTHON_CONVERSION::to_python(v[i])));
      }
      return t;
    }

    static
    cctbx::shared_storage<ValueType>&
    as_1d(cctbx::shared_storage_nd<ValueType, cctbx::dimension<3> >& v) {
      return v.as_1d();
    }

  };

  template <typename ValueType>
  void
  wrap_shared_storage(boost::python::module_builder& module,
                      const std::string& python_name,
                      const cctbx::type_holder<ValueType>)
  {
    boost::python::class_builder<
      cctbx::shared_storage<ValueType>,
      shared_storage_wrapper<ValueType> >
    py_ss(module, python_name.c_str());
    boost::python::export_converters(py_ss);

    boost::python::class_builder<
      cctbx::shared_storage_nd<ValueType, cctbx::dimension<3> > >
    py_ss_3d(module, (python_name + "_3d").c_str());
    boost::python::export_converters(py_ss_3d);

    py_ss_3d.declare_base(py_ss, boost::python::without_downcast);

    py_ss.def(boost::python::constructor<>());
    py_ss.def(boost::python::constructor<std::size_t>());
    py_ss.def(boost::python::constructor<
      const cctbx::shared_storage_handle_type&, std::size_t>());
    py_ss.def(boost::python::constructor<boost::python::tuple>());
    py_ss.def(&cctbx::shared_storage<ValueType>::size, "size");
    py_ss.def(&cctbx::shared_storage<ValueType>::size, "__len__");
    py_ss.def(&shared_storage_access<ValueType>::getitem, "__getitem__");
    py_ss.def(&shared_storage_access<ValueType>::setitem, "__setitem__");
    py_ss.def(&cctbx::shared_storage<ValueType>::deepcopy, "deepcopy");
    py_ss.def(&shared_storage_access<ValueType>::handle, "handle");
    py_ss.def(&shared_storage_access<ValueType>::fill, "fill");
    py_ss.def(&shared_storage_access<ValueType>::as_tuple, "as_tuple");

    py_ss_3d.def(boost::python::constructor<>());
    py_ss_3d.def(boost::python::constructor<const cctbx::dimension<3>&>());
    py_ss_3d.def(
      &cctbx::shared_storage_nd<ValueType, cctbx::dimension<3> >::dim, "dim");
    py_ss_3d.def(&shared_storage_access<ValueType>::as_1d, "as_1d");
  }

} // namespace cctbx

#endif // CCTBX_SHARED_STORAGE_BPL_H
