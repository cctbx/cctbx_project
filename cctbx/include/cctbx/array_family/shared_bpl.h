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

#include <cctbx/array_family/reductions.h>
#include <cctbx/array_family/operator_functors.h>
#include <cctbx/array_family/generic_array_functors.h>
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
    typedef ElementType e_t;
    typedef af::shared<ElementType> sh_t;
    typedef typename af::integer_to_float<ElementType>::float_type f_e_t;

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

    // This type is needed only to work around a Visual C++ 6 bug.
    typedef
    std::pair<
      boost::python::class_builder<
        shared<ElementType>,
        shared_wrapper<ElementType> >,
      boost::python::class_builder<
        shared_items<ElementType> >
    >
    sh_class_builders;

    static
    sh_class_builders
    plain(
      boost::python::module_builder& bpl_module,
      const std::string& python_name)
    {
      using namespace boost::python;

      class_builder<shared_items<ElementType> >
      py_shared_items(bpl_module, (python_name+"_items").c_str());

      class_builder<
        shared<ElementType>,
        shared_wrapper<ElementType> >
      py_shared(bpl_module, python_name.c_str());
      export_converters(py_shared);

      py_shared_items.def(constructor<>());
      py_shared_items.def(constructor<shared<ElementType> const&>());
      py_shared_items.def(&shared_items<ElementType>::size, "__len__");
      py_shared_items.def(&shared_items<ElementType>::getitem, "__getitem__");

      py_shared.def(constructor<>());
      py_shared.def(constructor<std::size_t>());
      py_shared.def(constructor<tuple>());
      py_shared.def(size, "size");
      py_shared.def(size, "__len__");
      py_shared.def(capacity, "capacity");
      py_shared.def(getitem, "__getitem__");
      py_shared.def(setitem, "__setitem__");
      py_shared.def(front, "front");
      py_shared.def(back, "back");
      py_shared.def(fill, "fill");
      py_shared.def(deep_copy, "deep_copy");
      py_shared.def(reserve, "reserve");
      py_shared.def(assign, "assign");
      py_shared.def(push_back, "push_back");
      py_shared.def(push_back, "append");
      py_shared.def(pop_back, "pop_back");
      py_shared.def(insert_i_x, "insert");
      py_shared.def(insert_i_n_x, "insert");
      py_shared.def(erase_i, "erase");
      py_shared.def(erase_i_j, "erase");
      py_shared.def(resize, "resize");
      py_shared.def(clear, "clear");
      py_shared.def(append, "append");
      py_shared.def(indices, "indices");
      py_shared.def(items, "items");

      return std::make_pair(py_shared, py_shared_items);
    }

    static
    shared<bool>
    invert_a(shared<bool> const& a)
    {
      shared<ElementType> result;
      result.reserve(a.size());
      for(std::size_t i=0;i<a.size();i++) result.append(!a[i]);
      return result;
    }

    static
    shared<bool>
    and_a_a(shared<bool> const& a1, shared<bool> const& a2)
    {
      if (a1.size() != a2.size()) throw_range_error();
      shared<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] && a2[i]);
      return result;
    }

    static
    shared<bool>
    or_a_a(shared<bool> const& a1, shared<bool> const& a2)
    {
      if (a1.size() != a2.size()) throw_range_error();
      shared<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] || a2[i]);
      return result;
    }

    static
    sh_class_builders
    logical(
      boost::python::module_builder& bpl_module,
      const std::string& python_name)
    {
      sh_class_builders class_blds = plain(bpl_module, python_name);
      class_blds.first.def(invert_a, "__invert__");
      class_blds.first.def(and_a_a, "__and__");
      class_blds.first.def(or_a_a, "__or__");
      return class_blds;
    }

    static
    shared<double>
    as_double(shared<ElementType> const& a)
    {
      shared<double> result;
      result.reserve(a.size());
      for(std::size_t i=0;i<a.size();i++) result.append(a[i]);
      return result;
    }

    static
    shared<ElementType>
    neg_a(shared<ElementType> const& a)
    {
      shared<ElementType> result;
      result.reserve(a.size());
      for(std::size_t i=0;i<a.size();i++) result.append(-a[i]);
      return result;
    }

    static
    shared<ElementType>
    add_a_a(shared<ElementType> const& a1, shared<ElementType> const& a2)
    {
      if (a1.size() != a2.size()) throw_range_error();
      shared<ElementType> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] + a2[i]);
      return result;
    }

    static
    shared<ElementType>
    sub_a_a(shared<ElementType> const& a1, shared<ElementType> const& a2)
    {
      if (a1.size() != a2.size()) throw_range_error();
      shared<ElementType> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] - a2[i]);
      return result;
    }

    static
    shared<ElementType>
    mul_a_a(shared<ElementType> const& a1, shared<ElementType> const& a2)
    {
      if (a1.size() != a2.size()) throw_range_error();
      shared<ElementType> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] * a2[i]);
      return result;
    }

    static
    shared<ElementType>
    div_a_a(shared<ElementType> const& a1, shared<ElementType> const& a2)
    {
      if (a1.size() != a2.size()) throw_range_error();
      shared<ElementType> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] / a2[i]);
      return result;
    }

    static
    shared<ElementType>
    mod_a_a(shared<ElementType> const& a1, shared<ElementType> const& a2)
    {
      if (a1.size() != a2.size()) throw_range_error();
      shared<ElementType> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] % a2[i]);
      return result;
    }

    static
    shared<ElementType>
    add_a_s(shared<ElementType> const& a1, ElementType const& a2)
    {
      shared<ElementType> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] + a2);
      return result;
    }

    static
    shared<ElementType>
    sub_a_s(shared<ElementType> const& a1, ElementType const& a2)
    {
      shared<ElementType> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] - a2);
      return result;
    }

    static
    shared<ElementType>
    mul_a_s(shared<ElementType> const& a1, ElementType const& a2)
    {
      shared<ElementType> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] * a2);
      return result;
    }

    static
    shared<ElementType>
    div_a_s(shared<ElementType> const& a1, ElementType const& a2)
    {
      shared<ElementType> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] / a2);
      return result;
    }

    static
    shared<ElementType>
    mod_a_s(shared<ElementType> const& a1, ElementType const& a2)
    {
      shared<ElementType> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] % a2);
      return result;
    }

    static
    shared<ElementType>
    iadd_a_s(shared<ElementType>& a1, ElementType const& a2)
    {
      for(std::size_t i=0;i<a1.size();i++) a1[i] += a2;
      return a1;
    }

    static
    shared<ElementType>
    isub_a_s(shared<ElementType>& a1, ElementType const& a2)
    {
      for(std::size_t i=0;i<a1.size();i++) a1[i] -= a2;
      return a1;
    }

    static
    shared<ElementType>
    imul_a_s(shared<ElementType>& a1, ElementType const& a2)
    {
      for(std::size_t i=0;i<a1.size();i++) a1[i] *= a2;
      return a1;
    }

    static
    shared<ElementType>
    idiv_a_s(shared<ElementType>& a1, ElementType const& a2)
    {
      for(std::size_t i=0;i<a1.size();i++) a1[i] /= a2;
      return a1;
    }

    static
    shared<ElementType>
    imod_a_s(shared<ElementType>& a1, ElementType const& a2)
    {
      for(std::size_t i=0;i<a1.size();i++) a1[i] %= a2;
      return a1;
    }

    static
    shared<bool>
    eq_a_a(shared<ElementType> const& a1, shared<ElementType> const& a2)
    {
      if (a1.size() != a2.size()) throw_range_error();
      shared<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] == a2[i]);
      return result;
    }

    static
    shared<bool>
    ne_a_a(shared<ElementType> const& a1, shared<ElementType> const& a2)
    {
      if (a1.size() != a2.size()) throw_range_error();
      shared<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] != a2[i]);
      return result;
    }

    static
    shared<bool>
    lt_a_a(shared<ElementType> const& a1, shared<ElementType> const& a2)
    {
      if (a1.size() != a2.size()) throw_range_error();
      shared<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] < a2[i]);
      return result;
    }

    static
    shared<bool>
    gt_a_a(shared<ElementType> const& a1, shared<ElementType> const& a2)
    {
      if (a1.size() != a2.size()) throw_range_error();
      shared<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] > a2[i]);
      return result;
    }

    static
    shared<bool>
    le_a_a(shared<ElementType> const& a1, shared<ElementType> const& a2)
    {
      if (a1.size() != a2.size()) throw_range_error();
      shared<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] <= a2[i]);
      return result;
    }

    static
    shared<bool>
    ge_a_a(shared<ElementType> const& a1, shared<ElementType> const& a2)
    {
      if (a1.size() != a2.size()) throw_range_error();
      shared<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] >= a2[i]);
      return result;
    }

    static
    shared<bool>
    eq_a_s(shared<ElementType> const& a1, ElementType const& a2)
    {
      shared<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] == a2);
      return result;
    }

    static
    shared<bool>
    ne_a_s(shared<ElementType> const& a1, ElementType const& a2)
    {
      shared<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] != a2);
      return result;
    }

    static
    shared<bool>
    lt_a_s(shared<ElementType> const& a1, ElementType const& a2)
    {
      shared<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] < a2);
      return result;
    }

    static
    shared<bool>
    gt_a_s(shared<ElementType> const& a1, ElementType const& a2)
    {
      shared<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] > a2);
      return result;
    }

    static
    shared<bool>
    le_a_s(shared<ElementType> const& a1, ElementType const& a2)
    {
      shared<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] <= a2);
      return result;
    }

    static
    shared<bool>
    ge_a_s(shared<ElementType> const& a1, ElementType const& a2)
    {
      shared<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] >= a2);
      return result;
    }

    static std::size_t
    max_index(sh_t const& a) { return af::max_index(a.const_ref()); }
    static std::size_t
    min_index(sh_t const& a) { return af::min_index(a.const_ref()); }
    static e_t max(sh_t const& a) { return af::max(a.const_ref()); }
    static e_t min(sh_t const& a) { return af::min(a.const_ref()); }
    static e_t sum(sh_t const& a) { return af::sum(a.const_ref()); }
    static e_t product(sh_t const& a) { return af::product(a.const_ref()); }
    static f_e_t mean(sh_t const& a) { return af::mean(a.const_ref()); }
    static f_e_t rms(sh_t const& a) { return af::rms(a.const_ref()); }

    static
    sh_class_builders
    numeric(
      boost::python::module_builder& bpl_module,
      const std::string& python_name)
    {
      sh_class_builders class_blds = plain(bpl_module, python_name);
      class_blds.first.def(as_double, "as_double");
      class_blds.first.def(neg_a, "__neg__");
      class_blds.first.def(add_a_a, "__add__");
      class_blds.first.def(sub_a_a, "__sub__");
      class_blds.first.def(mul_a_a, "__mul__");
      class_blds.first.def(div_a_a, "__div__");
      class_blds.first.def(add_a_s, "add");
      class_blds.first.def(sub_a_s, "sub");
      class_blds.first.def(mul_a_s, "mul");
      class_blds.first.def(div_a_s, "div");
      class_blds.first.def(iadd_a_s, "__iadd__");
      class_blds.first.def(isub_a_s, "__isub__");
      class_blds.first.def(imul_a_s, "__imul__");
      class_blds.first.def(idiv_a_s, "__idiv__");
      class_blds.first.def(eq_a_a, "__eq__");
      class_blds.first.def(ne_a_a, "__ne__");
      class_blds.first.def(lt_a_a, "__lt__");
      class_blds.first.def(gt_a_a, "__gt__");
      class_blds.first.def(le_a_a, "__le__");
      class_blds.first.def(ge_a_a, "__ge__");
      class_blds.first.def(eq_a_s, "__eq__");
      class_blds.first.def(ne_a_s, "__ne__");
      class_blds.first.def(lt_a_s, "__lt__");
      class_blds.first.def(gt_a_s, "__gt__");
      class_blds.first.def(le_a_s, "__le__");
      class_blds.first.def(ge_a_s, "__ge__");
      bpl_module.def(min_index, "min_index");
      bpl_module.def(max_index, "max_index");
      bpl_module.def(min, "min");
      bpl_module.def(max, "max");
      bpl_module.def(sum, "sum");
      bpl_module.def(product, "product");
      bpl_module.def(mean, "mean");
      bpl_module.def(rms, "rms");
      return class_blds;
    }

    static
    sh_class_builders
    integer(
      boost::python::module_builder& bpl_module,
      const std::string& python_name)
    {
      sh_class_builders class_blds = numeric(bpl_module, python_name);
      class_blds.first.def(mod_a_a, "__mod__");
      class_blds.first.def(mod_a_s, "mod");
      class_blds.first.def(imod_a_s, "__imod__");
      return class_blds;
    }
  };

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_SHARED_BPL_H
