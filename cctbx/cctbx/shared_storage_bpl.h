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

namespace boost {
  namespace python {

    // A wrapper is used to define additional constructors.
    template <typename T>
    struct shared_storage_wrapper : cctbx::shared_storage<T>
    {
      // Tell the compiler how to convert a base class object to
      // this wrapper object.
      shared_storage_wrapper(PyObject*,
                             const cctbx::shared_storage<T>& v)
        : cctbx::shared_storage<T>(v)
      {}
      shared_storage_wrapper(PyObject*)
        : cctbx::shared_storage<T>()
      {}
      shared_storage_wrapper(PyObject* self,
                             std::size_t n)
        : cctbx::shared_storage<T>(n)
      {}
      shared_storage_wrapper(PyObject* self,
                             const cctbx::shared_storage_handle_type& handle,
                             std::size_t n)
        : cctbx::shared_storage<T>(handle, n)
      {}
      shared_storage_wrapper(PyObject* self,
                             boost::python::tuple tuple)
        : cctbx::shared_storage<T>(tuple.size())
      {
        cctbx::shared_storage<T>::iterator v = this->begin();
        for (std::size_t i = 0; i < tuple.size(); i++)
          v[i] = BOOST_PYTHON_CONVERSION::from_python(tuple[i].get(),
                                              boost::python::type<T>());
      }

    };

    void raise_IndexError() {
      PyErr_SetString(PyExc_IndexError, "index out of range");
      throw boost::python::error_already_set();
    }

    template <typename T>
    struct shared_storage_access
    {
      static
      T
      getitem(const cctbx::shared_storage<T>& v,
              std::size_t key)
      {
        if (key >= v.size()) raise_IndexError();
        return v[key];
      }

      static
      void
      setitem(cctbx::shared_storage<T>& v,
              std::size_t key,
              const T &value)
      {
        if (key >= v.size()) raise_IndexError();
        v[key] = value;
      }

      static
      cctbx::shared_storage_handle_type&
      handle(cctbx::shared_storage<T>& v)
      {
        return v.handle();
      }

      static
      void fill(cctbx::shared_storage<T>& v, const T& value)
      {
        std::fill(v.begin(), v.end(), value);
      }

      // Convert to a regular Python tuple.
      static
      boost::python::tuple
      as_tuple(const cctbx::shared_storage<T>& v)
      {
        boost::python::tuple t(v.size());
        for (std::size_t i = 0; i < v.size(); i++) {
          t.set_item(i,
            boost::python::ref(BOOST_PYTHON_CONVERSION::to_python(v[i])));
        }
        return t;
      }

      static
      cctbx::shared_storage<T>&
      as_1d(cctbx::shared_storage_nd<T, cctbx::dimension<3> >& v) {
        return v.as_1d();
      }

    };

    template <typename T>
    void
    wrap_shared_storage(boost::python::module_builder& module,
                        const std::string& python_name,
                        const cctbx::type_holder<T>)
    {
      boost::python::class_builder<
        cctbx::shared_storage<T>,
        shared_storage_wrapper<T> >
      py_type(module, python_name.c_str());
      boost::python::export_converters(py_type);

      boost::python::class_builder<
        cctbx::shared_storage_nd<T, cctbx::dimension<3> > >
      py_type_3d(module, (python_name + "_3d").c_str());
      boost::python::export_converters(py_type_3d);

      py_type_3d.declare_base(py_type, boost::python::without_downcast);

      py_type.def(boost::python::constructor<>());
      py_type.def(boost::python::constructor<std::size_t>());
      py_type.def(boost::python::constructor<
        const cctbx::shared_storage_handle_type&, std::size_t>());
      py_type.def(boost::python::constructor<boost::python::tuple>());
      py_type.def(&cctbx::shared_storage<T>::size, "size");
      py_type.def(&cctbx::shared_storage<T>::size, "__len__");
      py_type.def(&shared_storage_access<T>::getitem, "__getitem__");
      py_type.def(&shared_storage_access<T>::setitem, "__setitem__");
      py_type.def(&cctbx::shared_storage<T>::deepcopy, "deepcopy");
      py_type.def(&shared_storage_access<T>::handle, "handle");
      py_type.def(&shared_storage_access<T>::fill, "fill");
      py_type.def(&shared_storage_access<T>::as_tuple, "as_tuple");

      py_type_3d.def(boost::python::constructor<>());
      py_type_3d.def(boost::python::constructor<const dimension<3>&>());
      py_type_3d.def(&cctbx::shared_storage_nd<T, dimension<3> >::dim, "dim");
      py_type_3d.def(&shared_storage_access<T>::as_1d, "as_1d");
    }

  } // namespace python
} // namespace boost

#endif // CCTBX_SHARED_STORAGE_BPL_H
