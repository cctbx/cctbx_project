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

#include <cctbx/constants.h>
#include <cctbx/math/utils.h>
#include <cctbx/array_family/reductions.h>
#include <cctbx/array_family/operator_functors.h>
#include <cctbx/array_family/generic_array_functors.h>
#include <cctbx/array_family/shared.h>

namespace cctbx { namespace af {

  inline void raise_IndexError() {
    PyErr_SetString(PyExc_IndexError, "index out of range");
    boost::python::throw_error_already_set();
  }

  inline void raise_incompatible_sizes() {
    PyErr_SetString(PyExc_RuntimeError, "incompatible array sizes");
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

  boost::python::ref shared_bool_getstate(shared<bool> const& a);
  void shared_bool_setstate(shared<bool>& a, boost::python::ref state);
  boost::python::ref shared_int_getstate(shared<int> const& a);
  void shared_int_setstate(shared<int>& a, boost::python::ref state);
  boost::python::ref shared_long_getstate(shared<long> const& a);
  void shared_long_setstate(shared<long>& a, boost::python::ref state);
  boost::python::ref shared_float_getstate(shared<float> const& a);
  void shared_float_setstate(shared<float>& a, boost::python::ref state);
  boost::python::ref shared_double_getstate(shared<double> const& a);
  void shared_double_setstate(shared<double>& a, boost::python::ref state);
  boost::python::ref shared_complex_double_getstate(
    shared<std::complex<double> > const& a);
  void shared_complex_double_setstate(
    shared<std::complex<double> >& a,
    boost::python::ref state);

  template <typename ElementType>
  struct shared_pickle
  {
    template <typename ClassBuilderType>
    static void def(ClassBuilderType& class_bldr) {}
  };

  template <>
  struct shared_pickle<bool>
  {
    template <typename ClassBuilderType>
    static void def(ClassBuilderType& class_bldr)
    {
      class_bldr.def(shared_bool_getstate, "__getstate__");
      class_bldr.def(shared_bool_setstate, "__setstate__");
    }
  };

  template <>
  struct shared_pickle<int>
  {
    template <typename ClassBuilderType>
    static void def(ClassBuilderType& class_bldr)
    {
      class_bldr.def(shared_int_getstate, "__getstate__");
      class_bldr.def(shared_int_setstate, "__setstate__");
    }
  };

  template <>
  struct shared_pickle<long>
  {
    template <typename ClassBuilderType>
    static void def(ClassBuilderType& class_bldr)
    {
      class_bldr.def(shared_long_getstate, "__getstate__");
      class_bldr.def(shared_long_setstate, "__setstate__");
    }
  };

  template <>
  struct shared_pickle<float>
  {
    template <typename ClassBuilderType>
    static void def(ClassBuilderType& class_bldr)
    {
      class_bldr.def(shared_float_getstate, "__getstate__");
      class_bldr.def(shared_float_setstate, "__setstate__");
    }
  };

  template <>
  struct shared_pickle<double>
  {
    template <typename ClassBuilderType>
    static void def(ClassBuilderType& class_bldr)
    {
      class_bldr.def(shared_double_getstate, "__getstate__");
      class_bldr.def(shared_double_setstate, "__setstate__");
    }
  };

  template <>
  struct shared_pickle<std::complex<double> >
  {
    template <typename ClassBuilderType>
    static void def(ClassBuilderType& class_bldr)
    {
      class_bldr.def(shared_complex_double_getstate, "__getstate__");
      class_bldr.def(shared_complex_double_setstate, "__setstate__");
    }
  };

  template <typename ElementType>
  struct shared_wrapper : shared<ElementType>
  {
    typedef ElementType e_t;
    typedef shared<ElementType> sh_t;

    typedef
    std::pair<
      boost::python::class_builder<
        sh_t,
        shared_wrapper<ElementType> >,
      boost::python::class_builder<
        shared_items<ElementType> >
    >
    sh_class_builders;

    // Tell the compiler how to convert a base class object to
    // this wrapper object.
    shared_wrapper(PyObject*, sh_t const& v)
      : sh_t(v)
    {}
    shared_wrapper(PyObject*)
      : sh_t()
    {}
    shared_wrapper(PyObject* self, std::size_t n)
      : sh_t(n)
    {}
    shared_wrapper(PyObject* self, std::size_t n, ElementType const& x)
      : sh_t(n, x)
    {}
    shared_wrapper(PyObject* self, boost::python::tuple tuple)
      : sh_t(tuple.size())
    {
      sh_t::iterator v = this->begin();
      for (std::size_t i = 0; i < tuple.size(); i++)
        v[i] = BOOST_PYTHON_CONVERSION::from_python(
          tuple[i].get(), boost::python::type<e_t const&>());
    }

    static
    std::size_t
    id(sh_t const& v) { return v.id(); }

    static
    std::size_t
    size(sh_t const& v) { return v.size(); }

    static
    std::size_t
    capacity(sh_t const& v) { return v.capacity(); }

    static
    e_t
    getitem(sh_t const& v, std::size_t i)
    {
      if (i >= v.size()) raise_IndexError();
      return v[i];
    }

    static
    void
    setitem(sh_t& v, std::size_t i, e_t const& x)
    {
      if (i >= v.size()) raise_IndexError();
      v[i] = x;
    }

    static
    e_t
    front(sh_t const& v)
    {
      if (v.size() == 0) raise_IndexError();
      return v.front();
    }

    static
    e_t
    back(sh_t const& v) {
      if (v.size() == 0) raise_IndexError();
      return v.back();
    }

    static
    void
    fill(sh_t& v, e_t const& x)
    {
      v.fill(x);
    }

    static
    void
    reserve(sh_t& v, std::size_t sz)
    {
      v.reserve(sz);
    }

    static
    sh_t
    deep_copy(sh_t const& v)
    {
      return v.deep_copy();
    }

    static
    void
    assign(sh_t& v, std::size_t sz, e_t const& x)
    {
      v.assign(sz, x);
    }

    static
    void push_back(sh_t& v, e_t const& x)
    {
      v.push_back(x);
    }

    static
    void
    pop_back(sh_t& v)
    {
      if (v.size() == 0) raise_IndexError();
      v.pop_back();
    }

    static
    void
    insert_i_x(sh_t& v, std::size_t i, e_t const& x)
    {
      if (i >= v.size()) raise_IndexError();
      v.insert(&v[i], x);
    }

    static
    void
    insert_i_n_x(sh_t& v, std::size_t i, std::size_t n, e_t const& x)
    {
      if (i >= v.size()) raise_IndexError();
      v.insert(&v[i], n, x);
    }

    static
    void
    erase_i(sh_t& v, std::size_t i)
    {
      if (i >= v.size()) raise_IndexError();
      v.erase(&v[i]);
    }

    static
    void erase_i_j(sh_t& v, std::size_t i, std::size_t j)
    {
      if (i >= v.size()) raise_IndexError();
      if (j >= v.size()) raise_IndexError();
      v.erase(&v[i], &v[j]);
    }

    static
    void
    resize(sh_t& v, std::size_t sz)
    {
      v.resize(sz);
    }

    static
    void
    clear(sh_t& v)
    {
      v.clear();
    }

    static
    void
    append(sh_t& v, sh_t& other)
    {
      v.insert(v.end(), other.begin(), other.end());
    }

    static
    boost::python::ref
    indices(sh_t const& v)
    {
      return boost::python::ref(PyRange_New(0, v.size(), 1, 1));
    }

    static
    shared_items<e_t>
    items(sh_t const& v)
    {
      return shared_items<e_t>(v);
    }

    static
    sh_t
    select(sh_t const& v, shared<bool> const& flags)
    {
      if (v.size() != flags.size()) raise_incompatible_sizes();
      std::size_t n = 0;
      std::size_t i;
      for(i=0;i<flags.size();i++) if (flags[i]) n++;
      sh_t result;
      result.reserve(n);
      for(i=0;i<flags.size();i++) if (flags[i]) result.push_back(v[i]);
      return result;
    }

    static
    shared<bool>
    invert_a(shared<bool> const& a)
    {
      sh_t result;
      result.reserve(a.size());
      for(std::size_t i=0;i<a.size();i++) result.append(!a[i]);
      return result;
    }

    static
    shared<bool>
    and_a_a(shared<bool> const& a1, shared<bool> const& a2)
    {
      if (a1.size() != a2.size()) raise_incompatible_sizes();
      shared<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] && a2[i]);
      return result;
    }

    static
    shared<bool>
    or_a_a(shared<bool> const& a1, shared<bool> const& a2)
    {
      if (a1.size() != a2.size()) raise_incompatible_sizes();
      shared<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] || a2[i]);
      return result;
    }

    static
    shared<bool>
    iand_a_a(shared<bool> a1, shared<bool> const& a2)
    {
      if (a1.size() != a2.size()) raise_incompatible_sizes();
      for(std::size_t i=0;i<a1.size();i++) if(!a2[i]) a1[i] = false;
      return a1;
    }

    static
    shared<bool>
    ior_a_a(shared<bool> a1, shared<bool> const& a2)
    {
      if (a1.size() != a2.size()) raise_incompatible_sizes();
      for(std::size_t i=0;i<a1.size();i++) if(a2[i]) a1[i] = true;
      return a1;
    }

    static
    shared<bool>
    iand_a_s(shared<bool> a1, bool a2)
    {
      if (!a2) std::fill(a1.begin(), a1.end(), false);
      return a1;
    }

    static
    shared<bool>
    ior_a_s(shared<bool> a1, bool a2)
    {
      if (a2) std::fill(a1.begin(), a1.end(), true);
      return a1;
    }

    static
    std::size_t
    count(shared<bool> const& a1, bool a2)
    {
      std::size_t result = 0;
      for(std::size_t i=0;i<a1.size();i++) if(a1[i] == a2) result++;
      return result;
    }

    static
    shared<double>
    as_double(sh_t const& a)
    {
      shared<double> result;
      result.reserve(a.size());
      for(std::size_t i=0;i<a.size();i++) result.append(a[i]);
      return result;
    }

    static
    sh_t
    neg_a(sh_t const& a)
    {
      sh_t result;
      result.reserve(a.size());
      for(std::size_t i=0;i<a.size();i++) result.append(-a[i]);
      return result;
    }

    static
    sh_t
    add_a_a(sh_t const& a1, sh_t const& a2)
    {
      if (a1.size() != a2.size()) raise_incompatible_sizes();
      sh_t result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] + a2[i]);
      return result;
    }

    static
    sh_t
    sub_a_a(sh_t const& a1, sh_t const& a2)
    {
      if (a1.size() != a2.size()) raise_incompatible_sizes();
      sh_t result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] - a2[i]);
      return result;
    }

    static
    sh_t
    mul_a_a(sh_t const& a1, sh_t const& a2)
    {
      if (a1.size() != a2.size()) raise_incompatible_sizes();
      sh_t result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] * a2[i]);
      return result;
    }

    static
    sh_t
    div_a_a(sh_t const& a1, sh_t const& a2)
    {
      if (a1.size() != a2.size()) raise_incompatible_sizes();
      sh_t result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] / a2[i]);
      return result;
    }

    static
    sh_t
    mod_a_a(sh_t const& a1, sh_t const& a2)
    {
      if (a1.size() != a2.size()) raise_incompatible_sizes();
      sh_t result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] % a2[i]);
      return result;
    }

    static
    sh_t
    add_a_s(sh_t const& a1, e_t const& a2)
    {
      sh_t result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] + a2);
      return result;
    }

    static
    sh_t
    sub_a_s(sh_t const& a1, e_t const& a2)
    {
      sh_t result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] - a2);
      return result;
    }

    static
    sh_t
    mul_a_s(sh_t const& a1, e_t const& a2)
    {
      sh_t result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] * a2);
      return result;
    }

    static
    sh_t
    div_a_s(sh_t const& a1, e_t const& a2)
    {
      sh_t result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] / a2);
      return result;
    }

    static
    sh_t
    mod_a_s(sh_t const& a1, e_t const& a2)
    {
      sh_t result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] % a2);
      return result;
    }

    static
    sh_t
    iadd_a_s(sh_t& a1, e_t const& a2)
    {
      for(std::size_t i=0;i<a1.size();i++) a1[i] += a2;
      return a1;
    }

    static
    sh_t
    isub_a_s(sh_t& a1, e_t const& a2)
    {
      for(std::size_t i=0;i<a1.size();i++) a1[i] -= a2;
      return a1;
    }

    static
    sh_t
    imul_a_s(sh_t& a1, e_t const& a2)
    {
      for(std::size_t i=0;i<a1.size();i++) a1[i] *= a2;
      return a1;
    }

    static
    sh_t
    idiv_a_s(sh_t& a1, e_t const& a2)
    {
      for(std::size_t i=0;i<a1.size();i++) a1[i] /= a2;
      return a1;
    }

    static
    sh_t
    imod_a_s(sh_t& a1, e_t const& a2)
    {
      for(std::size_t i=0;i<a1.size();i++) a1[i] %= a2;
      return a1;
    }

    static
    int
    cmp_a_a(sh_t const& a1, sh_t const& a2)
    {
      if (a1.id() == a2.id()) return 0;
      if (a1.size() < a2.size()) return -1;
      if (a1.size() > a2.size()) return  1;
      for(std::size_t i=0;i<a1.size();i++) {
        if (a1[i] < a2[i]) return -1;
        if (a1[i] > a2[i]) return  1;
      }
      return 0;
    }

    static
    int
    cmp_a_s(sh_t const& a1, e_t const& a2)
    {
      for(std::size_t i=0;i<a1.size();i++) {
        if (a1[i] < a2) return -1;
        if (a1[i] > a2) return  1;
      }
      return 0;
    }

    static
    shared<bool>
    eq_a_a(sh_t const& a1, sh_t const& a2)
    {
      if (a1.size() != a2.size()) raise_incompatible_sizes();
      shared<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] == a2[i]);
      return result;
    }

    static
    shared<bool>
    ne_a_a(sh_t const& a1, sh_t const& a2)
    {
      if (a1.size() != a2.size()) raise_incompatible_sizes();
      shared<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] != a2[i]);
      return result;
    }

    static
    shared<bool>
    lt_a_a(sh_t const& a1, sh_t const& a2)
    {
      if (a1.size() != a2.size()) raise_incompatible_sizes();
      shared<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] < a2[i]);
      return result;
    }

    static
    shared<bool>
    gt_a_a(sh_t const& a1, sh_t const& a2)
    {
      if (a1.size() != a2.size()) raise_incompatible_sizes();
      shared<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] > a2[i]);
      return result;
    }

    static
    shared<bool>
    le_a_a(sh_t const& a1, sh_t const& a2)
    {
      if (a1.size() != a2.size()) raise_incompatible_sizes();
      shared<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] <= a2[i]);
      return result;
    }

    static
    shared<bool>
    ge_a_a(sh_t const& a1, sh_t const& a2)
    {
      if (a1.size() != a2.size()) raise_incompatible_sizes();
      shared<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] >= a2[i]);
      return result;
    }

    static
    shared<bool>
    eq_a_s(sh_t const& a1, e_t const& a2)
    {
      shared<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] == a2);
      return result;
    }

    static
    shared<bool>
    ne_a_s(sh_t const& a1, e_t const& a2)
    {
      shared<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] != a2);
      return result;
    }

    static
    shared<bool>
    lt_a_s(sh_t const& a1, e_t const& a2)
    {
      shared<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] < a2);
      return result;
    }

    static
    shared<bool>
    gt_a_s(sh_t const& a1, e_t const& a2)
    {
      shared<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] > a2);
      return result;
    }

    static
    shared<bool>
    le_a_s(sh_t const& a1, e_t const& a2)
    {
      shared<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] <= a2);
      return result;
    }

    static
    shared<bool>
    ge_a_s(sh_t const& a1, e_t const& a2)
    {
      shared<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.append(a1[i] >= a2);
      return result;
    }

    static
    sh_t
    abs_a(sh_t const& a)
    {
      sh_t result;
      result.reserve(a.size());
      for(std::size_t i=0;i<a.size();i++) result.append(math::abs(a[i]));
      return result;
    }

    static
    sh_t
    fmod_a_s(sh_t const& a1, e_t const& a2)
    {
      sh_t result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) {
        result.append(std::fmod(a1[i], a2));
      }
      return result;
    }

    static
    sh_t
    pow_a_s(sh_t const& a, e_t const& exponent)
    {
      sh_t result;
      result.reserve(a.size());
      for(std::size_t i=0;i<a.size();i++) {
        result.append(std::pow(a[i], exponent));
      }
      return result;
    }

    static
    sh_t
    atan2_a_a(sh_t const& a1, sh_t const& a2)
    {
      if (a1.size() != a2.size()) raise_incompatible_sizes();
      sh_t result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) {
        result.append(std::atan2(a1[i], a2[i]));
      }
      return result;
    }

#define CCTBX_ARRAY_FAMILY_SHARED_BPL_CMATH_1ARG(func) \
    static \
    sh_t \
    func ## _a(sh_t const& a) \
    { \
      sh_t result; \
      result.reserve(a.size()); \
      for(std::size_t i=0;i<a.size();i++) { \
        result.append(std::func(a[i])); \
      } \
      return result; \
    }

CCTBX_ARRAY_FAMILY_SHARED_BPL_CMATH_1ARG(acos)
CCTBX_ARRAY_FAMILY_SHARED_BPL_CMATH_1ARG(cos)
CCTBX_ARRAY_FAMILY_SHARED_BPL_CMATH_1ARG(tan)
CCTBX_ARRAY_FAMILY_SHARED_BPL_CMATH_1ARG(asin)
CCTBX_ARRAY_FAMILY_SHARED_BPL_CMATH_1ARG(cosh)
CCTBX_ARRAY_FAMILY_SHARED_BPL_CMATH_1ARG(tanh)
CCTBX_ARRAY_FAMILY_SHARED_BPL_CMATH_1ARG(atan)
CCTBX_ARRAY_FAMILY_SHARED_BPL_CMATH_1ARG(exp)
CCTBX_ARRAY_FAMILY_SHARED_BPL_CMATH_1ARG(sin)
CCTBX_ARRAY_FAMILY_SHARED_BPL_CMATH_1ARG(fabs)
CCTBX_ARRAY_FAMILY_SHARED_BPL_CMATH_1ARG(log)
CCTBX_ARRAY_FAMILY_SHARED_BPL_CMATH_1ARG(sinh)
CCTBX_ARRAY_FAMILY_SHARED_BPL_CMATH_1ARG(ceil)
CCTBX_ARRAY_FAMILY_SHARED_BPL_CMATH_1ARG(floor)
CCTBX_ARRAY_FAMILY_SHARED_BPL_CMATH_1ARG(log10)
CCTBX_ARRAY_FAMILY_SHARED_BPL_CMATH_1ARG(sqrt)

    static std::size_t
    max_index_a(sh_t const& a) { return af::max_index(a.const_ref()); }
    static std::size_t
    min_index_a(sh_t const& a) { return af::min_index(a.const_ref()); }
    static e_t max_a(sh_t const& a) { return af::max(a.const_ref()); }
    static e_t min_a(sh_t const& a) { return af::min(a.const_ref()); }
    static e_t sum_a(sh_t const& a) { return af::sum(a.const_ref()); }
    static e_t sum_sq_a(sh_t const& a) { return af::sum_sq(a.const_ref()); }
    static e_t product_a(sh_t const& a) { return af::product(a.const_ref()); }
    static e_t mean_a(sh_t const& a) { return af::mean(a.const_ref()); }
    static e_t mean_sq_a(sh_t const& a) { return af::mean_sq(a.const_ref()); }

    static
    e_t
    mean_weighted_a_a(sh_t const& a1, sh_t const& a2)
    {
      return af::mean_weighted(a1.const_ref(), a2.const_ref());
    }

    static
    e_t
    mean_sq_weighted_a_a(sh_t const& a1, sh_t const& a2)
    {
      return af::mean_sq_weighted(a1.const_ref(), a2.const_ref());
    }

    static
    shared<double>
    abs_complex(shared<std::complex<double> > const& a)
    {
      shared<double> result;
      result.reserve(a.size());
      for(std::size_t i=0;i<a.size();i++) {
        result.push_back(std::abs(a[i]));
      }
      return result;
    }

    static
    shared<double>
    arg_complex_2(shared<std::complex<double> > const& a, bool deg)
    {
      shared<double> result;
      result.reserve(a.size());
      for(std::size_t i=0;i<a.size();i++) {
        result.push_back(std::arg(a[i]));
        if (deg) result[i] /= constants::pi_180;
      }
      return result;
    }

    static
    shared<double>
    arg_complex_1(shared<std::complex<double> > const& a)
    {
      return arg_complex_2(a, false);
    }

    static
    shared<double>
    norm_complex(shared<std::complex<double> > const& a)
    {
      shared<double> result;
      result.reserve(a.size());
      for(std::size_t i=0;i<a.size();i++) {
        result.push_back(std::norm(a[i]));
      }
      return result;
    }

    static
    shared<std::complex<double> >
    polar_complex_3(
      shared<double> const& rho,
      shared<double> const& theta,
      bool deg)
    {
      if (rho.size() != theta.size()) raise_incompatible_sizes();
      shared<std::complex<double> > result;
      result.reserve(rho.size());
      if (deg) {
        for(std::size_t i=0;i<rho.size();i++) {
          result.push_back(std::polar(rho[i], theta[i] * constants::pi_180));
        }
      }
      else {
        for(std::size_t i=0;i<rho.size();i++) {
          result.push_back(std::polar(rho[i], theta[i]));
        }
      }
      return result;
    }

    static
    shared<std::complex<double> >
    polar_complex_2(
      shared<double> const& rho,
      shared<double> const& theta)
    {
      return polar_complex_3(rho, theta, false);
    }

    static
    sh_class_builders
    plain(
      boost::python::module_builder& bpl_module,
      std::string const& python_name)
    {
      using namespace boost::python;

      class_builder<shared_items<ElementType> >
      py_shared_items(bpl_module, (python_name+"_items").c_str());

      class_builder<
        sh_t,
        shared_wrapper<e_t> >
      py_shared(bpl_module, python_name.c_str());
      export_converters(py_shared);

      py_shared_items.def(constructor<>());
      py_shared_items.def(constructor<sh_t const&>());
      py_shared_items.def(&shared_items<ElementType>::size, "__len__");
      py_shared_items.def(&shared_items<ElementType>::getitem, "__getitem__");

      py_shared.def(constructor<>());
      py_shared.def(constructor<std::size_t>());
      py_shared.def(constructor<std::size_t, ElementType const&>());
      py_shared.def(constructor<tuple>());
      py_shared.def(id, "id");
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
      py_shared.def(select, "select");

      shared_pickle<e_t>::def(py_shared);

      return std::make_pair(py_shared, py_shared_items);
    }

    static
    sh_class_builders
    eq_comparable(
      boost::python::module_builder& bpl_module,
      std::string const& python_name)
    {
      sh_class_builders class_blds = plain(bpl_module, python_name);
      class_blds.first.def(eq_a_a, "__eq__");
      class_blds.first.def(ne_a_a, "__ne__");
      class_blds.first.def(eq_a_s, "__eq__");
      class_blds.first.def(ne_a_s, "__ne__");
      return class_blds;
    }

    static
    sh_class_builders
    cmp_comparable(
      boost::python::module_builder& bpl_module,
      std::string const& python_name)
    {
      sh_class_builders class_blds = eq_comparable(bpl_module, python_name);
      class_blds.first.def(cmp_a_a, "__cmp__");
      class_blds.first.def(cmp_a_s, "__cmp__");
      return class_blds;
    }

    static
    sh_class_builders
    logical(
      boost::python::module_builder& bpl_module,
      std::string const& python_name)
    {
      sh_class_builders class_blds = cmp_comparable(bpl_module, python_name);
      class_blds.first.def(invert_a, "__invert__");
      class_blds.first.def(and_a_a, "__and__");
      class_blds.first.def(or_a_a, "__or__");
      class_blds.first.def(iand_a_a, "__iand__");
      class_blds.first.def(ior_a_a, "__ior__");
      class_blds.first.def(iand_a_s, "__iand__");
      class_blds.first.def(ior_a_s, "__ior__");
      class_blds.first.def(count, "count");
      return class_blds;
    }

    static
    sh_class_builders
    numeric_common(
      boost::python::module_builder& bpl_module,
      std::string const& python_name)
    {
      sh_class_builders class_blds = eq_comparable(bpl_module, python_name);
      class_blds.first.def(neg_a, "__neg__");
      class_blds.first.def(add_a_a, "__add__");
      class_blds.first.def(sub_a_a, "__sub__");
      class_blds.first.def(mul_a_a, "__mul__");
      class_blds.first.def(div_a_a, "__div__");
      class_blds.first.def(add_a_s, "add"); // XXX __add__ did not work
      class_blds.first.def(sub_a_s, "sub"); // XXX __sub__ did not work
      class_blds.first.def(mul_a_s, "mul"); // XXX __mul__ did not work
      class_blds.first.def(div_a_s, "div"); // XXX __div__ did not work
      class_blds.first.def(iadd_a_s, "__iadd__");
      class_blds.first.def(isub_a_s, "__isub__");
      class_blds.first.def(imul_a_s, "__imul__");
      class_blds.first.def(idiv_a_s, "__idiv__");
      bpl_module.def(sum_a, "sum");
      bpl_module.def(sum_sq_a, "sum_sq");
      bpl_module.def(product_a, "product");
      return class_blds;
    }

    static
    sh_class_builders
    numeric_no_pow(
      boost::python::module_builder& bpl_module,
      std::string const& python_name)
    {
      sh_class_builders class_blds = numeric_common(bpl_module, python_name);
      class_blds.first.def(as_double, "as_double");
      class_blds.first.def(cmp_a_a, "__cmp__");
      class_blds.first.def(cmp_a_s, "cmp"); // XXX __cmp__ did not work
      class_blds.first.def(lt_a_a, "__lt__");
      class_blds.first.def(gt_a_a, "__gt__");
      class_blds.first.def(le_a_a, "__le__");
      class_blds.first.def(ge_a_a, "__ge__");
      class_blds.first.def(lt_a_s, "__lt__");
      class_blds.first.def(gt_a_s, "__gt__");
      class_blds.first.def(le_a_s, "__le__");
      class_blds.first.def(ge_a_s, "__ge__");
      bpl_module.def(abs_a, "abs");
      bpl_module.def(min_index_a, "min_index");
      bpl_module.def(max_index_a, "max_index");
      bpl_module.def(min_a, "min");
      bpl_module.def(max_a, "max");
      return class_blds;
    }

    static
    sh_class_builders
    numeric(
      boost::python::module_builder& bpl_module,
      std::string const& python_name)
    {
      sh_class_builders class_blds = numeric_no_pow(bpl_module, python_name);
      bpl_module.def(pow_a_s, "pow");
      bpl_module.def(fmod_a_s, "fmod");
      bpl_module.def(atan2_a_a, "atan2");
      bpl_module.def(acos_a, "acos");
      bpl_module.def(cos_a, "cos");
      bpl_module.def(tan_a, "tan");
      bpl_module.def(asin_a, "asin");
      bpl_module.def(cosh_a, "cosh");
      bpl_module.def(tanh_a, "tanh");
      bpl_module.def(atan_a, "atan");
      bpl_module.def(exp_a, "exp");
      bpl_module.def(sin_a, "sin");
      bpl_module.def(fabs_a, "fabs");
      bpl_module.def(log_a, "log");
      bpl_module.def(sinh_a, "sinh");
      bpl_module.def(ceil_a, "ceil");
      bpl_module.def(floor_a, "floor");
      bpl_module.def(log10_a, "log10");
      bpl_module.def(sqrt_a, "sqrt");
      bpl_module.def(mean_a, "mean");
      bpl_module.def(mean_sq_a, "mean_sq");
      bpl_module.def(mean_weighted_a_a, "mean_weighted");
      bpl_module.def(mean_sq_weighted_a_a, "mean_sq_weighted");
      return class_blds;
    }

    static
    sh_class_builders
    integer(
      boost::python::module_builder& bpl_module,
      std::string const& python_name)
    {
      sh_class_builders class_blds = numeric_no_pow(bpl_module, python_name);
      class_blds.first.def(mod_a_a, "__mod__");
      class_blds.first.def(mod_a_s, "mod");
      class_blds.first.def(imod_a_s, "__imod__");
      return class_blds;
    }

    static
    sh_class_builders
    complex(
      boost::python::module_builder& bpl_module,
      std::string const& python_name)
    {
      sh_class_builders class_blds = numeric_common(bpl_module, python_name);
      bpl_module.def(abs_complex, "abs");
      bpl_module.def(arg_complex_2, "arg");
      bpl_module.def(arg_complex_1, "arg");
      bpl_module.def(norm_complex, "norm");
      bpl_module.def(polar_complex_3, "polar");
      bpl_module.def(polar_complex_2, "polar");
      return class_blds;
    }
  };

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_SHARED_BPL_H
