// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Mar 2002: modified copy of shared_bpl.h (rwgk)
     2002 Jan 14: Created, based on std_vector_bpl.h (rwgk)
 */

#ifndef CCTBX_ARRAY_FAMILY_SHARED_BPL_H
#define CCTBX_ARRAY_FAMILY_SHARED_BPL_H

#include <cctbx/array_family/shared.h>

namespace cctbx { namespace af {

  inline void raise_IndexError() {
    PyErr_SetString(PyExc_IndexError, "index out of range");
    boost::python::throw_error_already_set();
  }

  template <typename ElementType>
  struct shared_items
  {
    shared_items() {}

    shared_items(shared<ElementType> const& data)
    : data_(data)
    {}

    std::size_t size() const { return data_.size(); }

    boost::python::tuple
    getitem(std::size_t i) const
    {
      if (i >= data_.size()) raise_IndexError();
      return boost::python::tuple(
        boost::python::make_ref(i), boost::python::make_ref(data_[i]));
    }

    shared<ElementType> data_;
  };

  // A wrapper is used to define additional constructors.
  template <typename ElementType>
  struct shared_wrapper : shared<ElementType>
  {
    // Tell the compiler how to convert a base class object to
    // this wrapper object.
    shared_wrapper(PyObject*,
                   const shared<ElementType>& v)
      : shared<ElementType>(v)
    {}
    shared_wrapper(PyObject*)
      : shared<ElementType>()
    {}
    shared_wrapper(PyObject* self,
                   std::size_t n)
      : shared<ElementType>(n)
    {}
    shared_wrapper(PyObject* self,
                   boost::python::tuple tuple)
      : shared<ElementType>(tuple.size())
    {
      shared<ElementType>::iterator v = this->begin();
      for (std::size_t i = 0; i < tuple.size(); i++)
        v[i] = BOOST_PYTHON_CONVERSION::from_python(
          tuple[i].get(), boost::python::type<const ElementType&>());
    }

  };

  template <typename ElementType>
  struct shared_access
  {
    static
    std::size_t
    size(shared<ElementType>& v) { return v.size(); }

    static
    std::size_t
    capacity(shared<ElementType>& v) { return v.capacity(); }

    static
    ElementType
    getitem(const shared<ElementType>& v,
            std::size_t i)
    {
      if (i >= v.size()) raise_IndexError();
      return v[i];
    }

    static
    void
    setitem(shared<ElementType>& v,
            std::size_t i,
            const ElementType& x)
    {
      if (i >= v.size()) raise_IndexError();
      v[i] = x;
    }

    static
    ElementType
    front(shared<ElementType>& v) {
      if (v.size() == 0) raise_IndexError();
      return v.front();
    }

    static
    ElementType
    back(shared<ElementType>& v) {
      if (v.size() == 0) raise_IndexError();
      return v.back();
    }

    static
    void
    fill(shared<ElementType>& v,
         const ElementType& x) { v.fill(x); }

    static
    void
    reserve(shared<ElementType>& v,
            std::size_t sz) { v.reserve(sz); }

    static
    shared<ElementType>
    deep_copy(shared<ElementType>& v) { return v.deep_copy(); }

    static
    void
    assign(shared<ElementType>& v,
           std::size_t sz,
           const ElementType& x) { v.assign(sz, x); }

    static
    void
    push_back(shared<ElementType>& v,
              const ElementType& x) { v.push_back(x); }

    static
    void
    pop_back(shared<ElementType>& v) {
      if (v.size() == 0) raise_IndexError();
      v.pop_back();
    }

    static
    void
    insert_i_x(shared<ElementType>& v,
               std::size_t i,
               const ElementType& x) {
      if (i >= v.size()) raise_IndexError();
      v.insert(&v[i], x);
    }

    static
    void
    insert_i_n_x(shared<ElementType>& v,
                 std::size_t i,
                 std::size_t n,
                 const ElementType& x) {
      if (i >= v.size()) raise_IndexError();
      v.insert(&v[i], n, x);
    }

    static
    void
    erase_i(shared<ElementType>& v,
            std::size_t i) {
      if (i >= v.size()) raise_IndexError();
      v.erase(&v[i]);
    }

    static
    void
    erase_i_j(shared<ElementType>& v,
              std::size_t i,
              std::size_t j) {
      if (i >= v.size()) raise_IndexError();
      if (j >= v.size()) raise_IndexError();
      v.erase(&v[i], &v[j]);
    }

    static
    void
    resize(shared<ElementType>& v,
           std::size_t sz) { v.resize(sz); }

    static
    void
    clear(shared<ElementType>& v) { v.clear(); }

    static
    void
    append(shared<ElementType>& v, shared<ElementType>& other) {
      v.insert(v.end(), other.begin(), other.end());
    }

    static
    boost::python::ref
    indices(shared<ElementType> const& v) {
      return boost::python::ref(PyRange_New(0, v.size(), 1, 1));
    }

    static
    shared_items<ElementType>
    items(shared<ElementType> const& v) {
      return shared_items<ElementType>(v);
    }

  };

  template <typename ElementType>
  struct wrap_shared {
    // The return value is needed only to work around a Visual C++ 6 bug.
    static
    std::pair<
      boost::python::class_builder<
        shared<ElementType>,
        shared_wrapper<ElementType> >,
      boost::python::class_builder<
        shared_items<ElementType> >
    >
    run(
      boost::python::module_builder& module,
      const std::string& python_name)
    {
      using namespace boost::python;

      class_builder<shared_items<ElementType> >
      py_shared_items(module, (python_name+"_items").c_str());

      class_builder<
        shared<ElementType>,
        shared_wrapper<ElementType> >
      py_shared(module, python_name.c_str());
      export_converters(py_shared);

      py_shared_items.def(constructor<>());
      py_shared_items.def(constructor<shared<ElementType> const&>());
      py_shared_items.def(&shared_items<ElementType>::size, "__len__");
      py_shared_items.def(&shared_items<ElementType>::getitem, "__getitem__");

      py_shared.def(constructor<>());
      py_shared.def(constructor<std::size_t>());
      py_shared.def(constructor<tuple>());
      py_shared.def(&shared_access<ElementType>::size, "size");
      py_shared.def(&shared_access<ElementType>::size, "__len__");
      py_shared.def(&shared_access<ElementType>::capacity, "capacity");
      py_shared.def(&shared_access<ElementType>::getitem, "__getitem__");
      py_shared.def(&shared_access<ElementType>::setitem, "__setitem__");
      py_shared.def(&shared_access<ElementType>::front, "front");
      py_shared.def(&shared_access<ElementType>::back, "back");
      py_shared.def(&shared_access<ElementType>::fill, "fill");
      py_shared.def(&shared_access<ElementType>::deep_copy, "deep_copy");
      py_shared.def(&shared_access<ElementType>::reserve, "reserve");
      py_shared.def(&shared_access<ElementType>::assign, "assign");
      py_shared.def(&shared_access<ElementType>::push_back, "push_back");
      py_shared.def(&shared_access<ElementType>::push_back, "append");
      py_shared.def(&shared_access<ElementType>::pop_back, "pop_back");
      py_shared.def(&shared_access<ElementType>::insert_i_x, "insert");
      py_shared.def(&shared_access<ElementType>::insert_i_n_x, "insert");
      py_shared.def(&shared_access<ElementType>::erase_i, "erase");
      py_shared.def(&shared_access<ElementType>::erase_i_j, "erase");
      py_shared.def(&shared_access<ElementType>::resize, "resize");
      py_shared.def(&shared_access<ElementType>::clear, "clear");
      py_shared.def(&shared_access<ElementType>::append, "append");
      py_shared.def(&shared_access<ElementType>::indices, "indices");
      py_shared.def(&shared_access<ElementType>::items, "items");

      return std::make_pair(py_shared, py_shared_items);
    }
  };

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_SHARED_BPL_H
