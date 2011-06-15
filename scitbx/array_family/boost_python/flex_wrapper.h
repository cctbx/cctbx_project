#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_FLEX_WRAPPER_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_FLEX_WRAPPER_H

#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_non_const_reference.hpp>
#include <boost/python/slice.hpp>
#include <boost_adaptbx/type_id_eq.h>
#include <scitbx/constants.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/versa_reductions.h>
#include <scitbx/array_family/versa_algebra.h>
#include <scitbx/boost_python/utils.h>
#include <scitbx/boost_python/slice.h>
#include <scitbx/boost_python/container_conversions.h>
#include <scitbx/array_family/boost_python/shared_flex_conversions.h>
#include <scitbx/array_family/boost_python/ref_flex_conversions.h>
#include <scitbx/array_family/boost_python/passing_flex_by_reference.h>
#include <scitbx/array_family/boost_python/selections_wrapper.h>
#include <scitbx/misc/positive_getitem_index.h>
#include <scitbx/array_family/slice.h>

namespace scitbx { namespace af { namespace boost_python {

  using scitbx::positive_getitem_index;

  template <typename ElementType>
  struct flex_default_element
  {
    static ElementType
    get() { return ElementType(); }
  };

  template <typename ElementType,
            typename GetitemReturnValuePolicy
              = boost::python::return_value_policy<
                  boost::python::copy_non_const_reference> >
  struct flex_wrapper : versa<ElementType, flex_grid<> >
  {
    typedef ElementType e_t;
    typedef versa<ElementType, flex_grid<> > f_t;
    typedef typename f_t::base_array_type base_array_type;

    // Tell the compiler how to convert a base class object to
    // this wrapper object.
    flex_wrapper(PyObject*, f_t const& v)
      : f_t(v)
    {}

    flex_wrapper(PyObject*)
      : f_t(flex_grid<>(0), flex_default_element<ElementType>::get())
    {}

    flex_wrapper(PyObject*, flex_grid<> const& fg)
      : f_t(fg, flex_default_element<ElementType>::get())
    {}

    flex_wrapper(PyObject*, flex_grid<> const& fg, ElementType const& x)
      : f_t(fg, x)
    {}

    flex_wrapper(PyObject*, std::size_t n)
      : f_t(flex_grid<>(n), flex_default_element<ElementType>::get())
    {}

    flex_wrapper(PyObject*, std::size_t n, ElementType const& x)
      : f_t(flex_grid<>(n), x)
    {}

    flex_wrapper(PyObject*, shared_plain<ElementType> const& a)
      : f_t(a, flex_grid<>(a.size()))
    {}

    static std::size_t
    element_size() { return sizeof(ElementType); }

    static flex_grid<>
    accessor(f_t const& a) { return a.accessor(); }

    static std::size_t
    nd(f_t const& a) { return a.accessor().nd(); }

    static bool
    is_0_based(f_t const& a) { return a.accessor().is_0_based(); }

    static flex_grid_default_index_type
    origin(f_t const& a) { return a.accessor().origin(); }

    static flex_grid_default_index_type
    all(f_t const& a) { return a.accessor().all(); }

    static flex_grid_default_index_type
    last_0(f_t const& a) { return a.accessor().last(); }

    static flex_grid_default_index_type
    last_1(f_t const& a, bool open_range)
    {
      return a.accessor().last(open_range);
    }

    static bool
    is_padded(f_t const& a) { return a.accessor().is_padded(); }

    static flex_grid_default_index_type
    focus_0(f_t const& a) { return a.accessor().focus(); }

    static flex_grid_default_index_type
    focus_1(f_t const& a, bool open_range)
    {
      return a.accessor().focus(open_range);
    }

    static std::size_t
    focus_size_1d(f_t const& a) { return a.accessor().focus_size_1d(); }

    static bool
    is_trivial_1d(f_t const& a) { return a.accessor().is_trivial_1d(); }

    static f_t
    shift_origin(f_t const& a)
    {
      return f_t(a, a.accessor().shift_origin());
    }

    static std::size_t
    id(f_t const& a) { return a.id(); }

    static std::size_t
    size(f_t const& a) { return a.size(); }

    static std::size_t
    capacity(f_t const& a) { return a.capacity(); }

    static e_t&
    getitem_1d(f_t& a, long i)
    {
      if (!a.check_shared_size()) raise_shared_size_mismatch();
      std::size_t j = positive_getitem_index(i, a.size());
      return a[j];
    }

    static f_t
    getitem_1d_slice(f_t const& a, boost::python::slice const& slice)
    {
      if (!a.check_shared_size()) raise_shared_size_mismatch();
      scitbx::boost_python::adapted_slice a_sl(slice, a.size());
      af::shared_plain<e_t> result(af::reserve(a_sl.size));
      for(long i=a_sl.start;i!=a_sl.stop;i+=a_sl.step) {
        result.push_back(a[i]);
      }
      return f_t(result, flex_grid<>(result.size()));
    }

    static e_t&
    getitem_flex_grid(f_t& a, flex_grid_default_index_type const& i)
    {
      if (!a.check_shared_size()) raise_shared_size_mismatch();
      if (!a.accessor().is_valid_index(i)) {
        scitbx::boost_python::raise_index_error();
      }
      return a(i);
    }

    static void
    setitem_1d(f_t& a, long i, e_t const& x)
    {
      if (!a.check_shared_size()) raise_shared_size_mismatch();
      std::size_t j = positive_getitem_index(i, a.size());
      a[j] = x;
    }

    static void
    setitem_flex_grid(
      f_t& a, flex_grid_default_index_type const& i, e_t const& x)
    {
      if (!a.check_shared_size()) raise_shared_size_mismatch();
      if (!a.accessor().is_valid_index(i)) {
        scitbx::boost_python::raise_index_error();
      }
      a(i) = x;
    }

    static void
    delitem_1d(f_t& a, long i)
    {
      base_array_type b = flex_as_base_array(a);
      std::size_t j = positive_getitem_index(i, b.size());
      b.erase(b.begin()+j);
      a.resize(
        flex_grid<>(b.size()), flex_default_element<ElementType>::get());
    }

    static void
    delitem_1d_slice(f_t& a, boost::python::slice const& slice)
    {
      base_array_type b = flex_as_base_array(a);
      scitbx::boost_python::adapted_slice a_sl(slice, b.size());
      SCITBX_ASSERT(a_sl.step == 1); // not implemented
      e_t* bb = b.begin();
      b.erase(bb+a_sl.start, bb+a_sl.stop);
      a.resize(
        flex_grid<>(b.size()), flex_default_element<ElementType>::get());
    }

    static e_t&
    front(f_t& a)
    {
      if (!a.check_shared_size()) raise_shared_size_mismatch();
      if (a.size() == 0) scitbx::boost_python::raise_index_error();
      return a.front();
    }

    static e_t&
    back(f_t& a)
    {
      if (!a.check_shared_size()) raise_shared_size_mismatch();
      if (a.size() == 0) scitbx::boost_python::raise_index_error();
      return a.back();
    }

    static void
    fill(f_t& a, e_t const& x)
    {
      if (!a.check_shared_size()) raise_shared_size_mismatch();
      a.fill(x);
    }

    static void
    reserve(f_t& a, std::size_t sz)
    {
      if (!a.check_shared_size()) raise_shared_size_mismatch();
      a.reserve(sz);
    }

    static f_t
    deep_copy(f_t const& a)
    {
      if (!a.check_shared_size()) raise_shared_size_mismatch();
      return a.deep_copy();
    }

    static f_t
    shallow_copy(f_t const& a)
    {
      if (!a.check_shared_size()) raise_shared_size_mismatch();
      return a;
    }

    static f_t
    as_1d(f_t const& a)
    {
      if (!a.check_shared_size()) raise_shared_size_mismatch();
      SCITBX_ASSERT(!a.accessor().is_padded());
      return f_t(a, flex_grid<>(a.size()));
    }

    static void
    assign(f_t& a, std::size_t sz, e_t const& x)
    {
      base_array_type b = flex_as_base_array(a);
      b.assign(sz, x);
      a.resize(
        flex_grid<>(b.size()), flex_default_element<ElementType>::get());
    }

    static void
    append(f_t& a, e_t const& x)
    {
      base_array_type b = flex_as_base_array(a);
      b.push_back(x);
      a.resize(
        flex_grid<>(b.size()), flex_default_element<ElementType>::get());
    }

    static void
    pop_back(f_t& a)
    {
      base_array_type b = flex_as_base_array(a);
      if (b.size() == 0) scitbx::boost_python::raise_index_error();
      b.pop_back();
      a.resize(
        flex_grid<>(b.size()), flex_default_element<ElementType>::get());
    }

    static void
    insert_i_x(f_t& a, long i, e_t const& x)
    {
      base_array_type b = flex_as_base_array(a);
      std::size_t j = positive_getitem_index(i, b.size(), true);
      b.insert(b.begin()+j, x);
      a.resize(
        flex_grid<>(b.size()), flex_default_element<ElementType>::get());
    }

    static void
    insert_i_n_x(f_t& a, long i, std::size_t n, e_t const& x)
    {
      base_array_type b = flex_as_base_array(a);
      std::size_t j = positive_getitem_index(i, b.size(), true);
      b.insert(b.begin()+j, n, x);
      a.resize(
        flex_grid<>(b.size()), flex_default_element<ElementType>::get());
    }

    static void
    resize_1d_1(f_t& a, std::size_t sz)
    {
      base_array_type b = flex_as_base_array(a);
      b.resize(sz, flex_default_element<ElementType>::get());
      a.resize(
        flex_grid<>(b.size()), flex_default_element<ElementType>::get());
    }

    static void
    resize_1d_2(f_t& a, std::size_t sz, e_t const& x)
    {
      base_array_type b = flex_as_base_array(a);
      b.resize(sz, x);
      a.resize(
        flex_grid<>(b.size()), flex_default_element<ElementType>::get());
    }

    static void
    resize_flex_grid_1(f_t& a, flex_grid<> const& grid)
    {
      a.resize(grid, flex_default_element<ElementType>::get());
    }

    static void
    resize_flex_grid_2(f_t& a, flex_grid<> const& grid, e_t const& x)
    {
      a.resize(grid, x);
    }

    static void
    reshape(f_t& a, flex_grid<> const& grid)
    {
      SCITBX_ASSERT(grid.size_1d() == a.size());
      a.resize(grid, flex_default_element<ElementType>::get());
    }

    static void
    clear(f_t& a)
    {
      base_array_type b = flex_as_base_array(a);
      b.clear();
      a.resize(
        flex_grid<>(b.size()), flex_default_element<ElementType>::get());
    }

    static void
    extend(f_t& a, f_t const& other)
    {
      base_array_type b = flex_as_base_array(a);
      assert_0_based_1d(other.accessor());
      b.insert(b.end(), other.begin(), other.end());
      a.resize(
        flex_grid<>(b.size()), flex_default_element<ElementType>::get());
    }

    static shared<e_t>
    concatenate(
      af::const_ref<e_t> const& self,
      af::const_ref<e_t> const& other)
    {
      shared<e_t> result((af::reserve(self.size()+other.size())));
      result.insert(result.end(), self.begin(), self.end());
      result.insert(result.end(), other.begin(), other.end());
      return result;
    }

    static shared<e_t>
    reversed(
      af::const_ref<e_t> const& self)
    {
      shared<e_t> result((af::reserve(self.size())));
      for(std::size_t i=self.size();i>0;) {
        i--;
        result.push_back(self[i]);
      }
      return result;
    }

    static boost::python::object
    set_selected_bool_a(
      boost::python::object flex_object,
      af::const_ref<bool> const& flags,
      af::const_ref<e_t> const& new_values)
    {
      boost::python::extract<af::ref<e_t> > a_proxy(flex_object);
      af::ref<e_t> a = a_proxy();
      SCITBX_ASSERT(a.size() == flags.size());
      if (new_values.size() == flags.size()) {
        e_t* ai = a.begin();
        const bool* fi = flags.begin();
        const e_t* ni = new_values.begin();
        const e_t* ne = new_values.end();
        while (ni != ne) {
          if (*fi++) *ai = *ni;
          ai++;
          ni++;
        }
      }
      else {
        std::size_t i_new_value = 0;
        for(std::size_t i=0;i<flags.size();i++) {
          if (flags[i]) {
            SCITBX_ASSERT(i_new_value < new_values.size());
            a[i] = new_values[i_new_value];
            i_new_value++;
          }
        }
        SCITBX_ASSERT(i_new_value == new_values.size());
      }
      return flex_object;
    }

    static boost::python::object
    set_selected_bool_s(
      boost::python::object flex_object,
      af::const_ref<bool, flex_grid<> > const& flags,
      e_t const& new_value)
    {
      boost::python::extract<af::ref<e_t, flex_grid<> > > a_proxy(flex_object);
      af::ref<e_t, flex_grid<> > a = a_proxy();
      SCITBX_ASSERT(a.accessor() == flags.accessor());
      for(std::size_t i=0;i<flags.size();i++) {
        if (flags[i]) a[i] = new_value;
      }
      return flex_object;
    }

    template <typename UnsignedType>
    static boost::python::object
    set_selected_unsigned_a(
      boost::python::object const& flex_object,
      af::const_ref<UnsignedType> const& indices,
      af::const_ref<e_t> const& new_values)
    {
      boost::python::extract<af::ref<e_t> > a_proxy(flex_object);
      af::ref<e_t> a = a_proxy();
      SCITBX_ASSERT(indices.size() == new_values.size());
      for(std::size_t i=0;i<indices.size();i++) {
        SCITBX_ASSERT(indices[i] < a.size());
        a[indices[i]] = new_values[i];
      }
      return flex_object;
    }

    template <typename UnsignedType>
    static boost::python::object
    set_selected_unsigned_s(
      boost::python::object const& flex_object,
      af::const_ref<UnsignedType> const& indices,
      e_t const& new_value)
    {
      boost::python::extract<af::ref<e_t> > a_proxy(flex_object);
      af::ref<e_t> a = a_proxy();
      for(std::size_t i=0;i<indices.size();i++) {
        SCITBX_ASSERT(indices[i] < a.size());
        a[indices[i]] = new_value;
      }
      return flex_object;
    }

    template <typename UnsignedType>
    static boost::python::object
    copy_selected_unsigned_a(
      boost::python::object const& flex_object,
      af::const_ref<UnsignedType> const& indices,
      af::const_ref<e_t> const& new_values)
    {
      boost::python::extract<af::ref<e_t> > a_proxy(flex_object);
      af::ref<e_t> a = a_proxy();
      SCITBX_ASSERT(a.size() == new_values.size());
      for(std::size_t i=0;i<indices.size();i++) {
        SCITBX_ASSERT(indices[i] < a.size());
        UnsignedType ii = indices[i];
        a[ii] = new_values[ii];
      }
      return flex_object;
    }

    static std::size_t
    count(f_t const& a1, e_t const& a2)
    {
      std::size_t result = 0;
      for(std::size_t i=0;i<a1.size();i++) if(a1[i] == a2) result++;
      return result;
    }

    static flex_double
    as_double(f_t const& a)
    {
      shared_plain<double> result(a.size(), init_functor_null<double>());
      for(std::size_t i=0;i<a.size();i++) result[i] = a[i];
      return flex_double(result, a.accessor());
    }

    static boost::optional<std::size_t>
    first_index_a_s(f_t const& a, e_t x) {
      return first_index(a, std::bind2nd(std::equal_to<e_t>(), x));
    }

    static boost::optional<std::size_t>
    last_index_a_s(f_t const& a, e_t x) {
      return last_index(a, std::bind2nd(std::equal_to<e_t>(), x));
    }

    static f_t neg_a(f_t const& a) { return -a; }
    static f_t add_a_a(f_t const& a1, f_t const& a2) { return a1 + a2; }
    static f_t sub_a_a(f_t const& a1, f_t const& a2) { return a1 - a2; }
    static f_t mul_a_a(f_t const& a1, f_t const& a2) { return a1 * a2; }
    static f_t div_a_a(f_t const& a1, f_t const& a2) { return a1 / a2; }
    static f_t mod_a_a(f_t const& a1, f_t const& a2) { return a1 % a2; }
    static f_t add_a_s(f_t const& a1, e_t const& a2) { return a1 + a2; }
    static f_t sub_a_s(f_t const& a1, e_t const& a2) { return a1 - a2; }
    static f_t rsub_a_s(f_t const& a2, e_t const& a1) { return a1 - a2; }
    static f_t mul_a_s(f_t const& a1, e_t const& a2) { return a1 * a2; }
    static f_t div_a_s(f_t const& a1, e_t const& a2) { return a1 / a2; }
    static f_t rdiv_a_s(f_t const& a2, e_t const& a1) { return a1 / a2; }
    static f_t mod_a_s(f_t const& a1, e_t const& a2) { return a1 % a2; }
    static f_t rmod_a_s(f_t const& a2, e_t const& a1) { return a1 % a2; }
    static f_t iadd_a_a(f_t& a1, f_t const& a2) { return a1 += a2; }
    static f_t isub_a_a(f_t& a1, f_t const& a2) { return a1 -= a2; }
    static f_t imul_a_a(f_t& a1, f_t const& a2) { return a1 *= a2; }
    static f_t idiv_a_a(f_t& a1, f_t const& a2) { return a1 /= a2; }
    static f_t imod_a_a(f_t& a1, f_t const& a2) { return a1 %= a2; }
    static f_t iadd_a_s(f_t& a1, e_t const& a2) { return a1 += a2; }
    static f_t isub_a_s(f_t& a1, e_t const& a2) { return a1 -= a2; }
    static f_t imul_a_s(f_t& a1, e_t const& a2) { return a1 *= a2; }
    static f_t idiv_a_s(f_t& a1, e_t const& a2) { return a1 /= a2; }
    static f_t imod_a_s(f_t& a1, e_t const& a2) { return a1 %= a2; }

    static bool
    all_eq_a_a(f_t const& a1, f_t const& a2) { return a1.all_eq(a2); }
    static bool
    all_ne_a_a(f_t const& a1, f_t const& a2) { return a1.all_ne(a2); }
    static bool
    all_lt_a_a(f_t const& a1, f_t const& a2) { return a1.all_lt(a2); }
    static bool
    all_gt_a_a(f_t const& a1, f_t const& a2) { return a1.all_gt(a2); }
    static bool
    all_le_a_a(f_t const& a1, f_t const& a2) { return a1.all_le(a2); }
    static bool
    all_ge_a_a(f_t const& a1, f_t const& a2) { return a1.all_ge(a2); }

    static bool
    all_eq_a_s(f_t const& a1, e_t const& a2) { return a1.all_eq(a2); }
    static bool
    all_ne_a_s(f_t const& a1, e_t const& a2) { return a1.all_ne(a2); }
    static bool
    all_lt_a_s(f_t const& a1, e_t const& a2) { return a1.all_lt(a2); }
    static bool
    all_gt_a_s(f_t const& a1, e_t const& a2) { return a1.all_gt(a2); }
    static bool
    all_le_a_s(f_t const& a1, e_t const& a2) { return a1.all_le(a2); }
    static bool
    all_ge_a_s(f_t const& a1, e_t const& a2) { return a1.all_ge(a2); }

    static int
    order_a_a(f_t const& a1, f_t const& a2) { return order(a1, a2); }

    static flex_bool eq_a_a(f_t const& a1, f_t const& a2) { return a1 == a2; }
#if defined(__GNUC__) && __GNUC__ == 2 && __GNUC_MINOR__ == 96
    static flex_bool ne_a_a(f_t const& a1, f_t const& a2)
    {
      if (a1.size() != a2.size()) throw_range_error();
      shared_plain<bool> result(a1.size(), init_functor_null<bool>());
      for(std::size_t i=0;i<a1.size();i++) result[i] = a1[i] != a2[1];
      return versa<bool, flex_grid<> >(result, a1.accessor());
    }
#else
    static flex_bool ne_a_a(f_t const& a1, f_t const& a2) { return a1 != a2; }
#endif
    static flex_bool lt_a_a(f_t const& a1, f_t const& a2) { return a1 < a2; }
    static flex_bool gt_a_a(f_t const& a1, f_t const& a2) { return a1 > a2; }
    static flex_bool le_a_a(f_t const& a1, f_t const& a2) { return a1 <= a2; }
    static flex_bool ge_a_a(f_t const& a1, f_t const& a2) { return a1 >= a2; }
    static flex_bool eq_a_s(f_t const& a1, e_t const& a2) { return a1 == a2; }
    static flex_bool ne_a_s(f_t const& a1, e_t const& a2) { return a1 != a2; }
    static flex_bool lt_a_s(f_t const& a1, e_t const& a2) { return a1 < a2; }
    static flex_bool gt_a_s(f_t const& a1, e_t const& a2) { return a1 > a2; }
    static flex_bool le_a_s(f_t const& a1, e_t const& a2) { return a1 <= a2; }
    static flex_bool ge_a_s(f_t const& a1, e_t const& a2) { return a1 >= a2; }

    static f_t abs_a(f_t const& a) { return absolute(a); }
    static f_t pow2_a(f_t const& a) { return pow2(a); }

    static f_t
    fmod_a_s(f_t const& a1, e_t const& a2) { return fmod(a1, a2); }

    static f_t
    fmod_positive_a_s(f_t const& a1, e_t const& a2)
    {
      return fmod_positive(a1, a2);
    }

    static f_t
    pow_a_s(f_t const& a, e_t const& exponent) { return pow(a, exponent); }

    static f_t
    atan2_a_a(f_t const& a1, f_t const& a2) { return atan2(a1, a2); }

    static f_t acos_a(f_t const& a) { return acos(a); }
    static f_t cos_a(f_t const& a) { return cos(a); }
    static f_t tan_a(f_t const& a) { return tan(a); }
    static f_t asin_a(f_t const& a) { return asin(a); }
    static f_t cosh_a(f_t const& a) { return cosh(a); }
    static f_t tanh_a(f_t const& a) { return tanh(a); }
    static f_t atan_a(f_t const& a) { return atan(a); }
    static f_t exp_a(f_t const& a) { return exp(a); }
    static f_t sin_a(f_t const& a) { return sin(a); }
    static f_t fabs_a(f_t const& a) { return fabs(a); }
    static f_t log_a(f_t const& a) { return log(a); }
    static f_t sinh_a(f_t const& a) { return sinh(a); }
    static f_t ceil_a(f_t const& a) { return ceil(a); }
    static f_t floor_a(f_t const& a) { return floor(a); }
    static f_t log10_a(f_t const& a) { return log10(a); }
    static f_t sqrt_a(f_t const& a) { return sqrt(a); }

    static std::size_t
    max_index_a(f_t const& a) { return max_index(a); }
    static std::size_t
    min_index_a(f_t const& a) { return min_index(a); }
    static e_t max_a(f_t const& a) { return max(a); }
    static e_t min_a(f_t const& a) { return min(a); }
    static e_t max_absolute_a(f_t const& a) { return max_absolute(a); }
    static e_t sum_a(f_t const& a) { return sum(a); }
    static e_t sum_sq_a(f_t const& a) { return sum_sq(a); }
    static e_t product_a(f_t const& a) { return product(a); }
    static e_t norm_a(f_t const& a) { return norm(a.const_ref()); }
    static e_t mean_a(f_t const& a) { return mean(a); }
    static e_t mean_sq_a(f_t const& a) { return mean_sq(a); }

    static e_t
    mean_weighted_a_a(f_t const& a1, f_t const& a2)
    {
      return mean_weighted(a1, a2);
    }

    static e_t
    mean_sq_weighted_a_a(f_t const& a1, f_t const& a2)
    {
      return mean_sq_weighted(a1, a2);
    }

    typedef boost::python::class_<f_t, flex_wrapper<ElementType> > class_f_t;

    static class_f_t
    plain(std::string const& python_name)
    {
      using namespace boost::python;

      scitbx::boost_python::container_conversions::from_python_sequence<
        shared_plain<ElementType>,
        scitbx::boost_python::container_conversions
          ::variable_capacity_policy>();
      scitbx::boost_python::container_conversions::from_python_sequence<
        shared<ElementType>,
        scitbx::boost_python::container_conversions
          ::variable_capacity_policy>();
      shared_flex_conversions<ElementType>();
      ref_flex_conversions<ElementType>();
      flex_1d_from_flex<ElementType>();

      class_f_t result(python_name.c_str());
      result
        .def(init<flex_grid<> const&, optional<ElementType const&> >())
        .def(init<std::size_t, optional<ElementType const&> >())
        .def(init<shared_plain<ElementType> const&>())
        .def("element_size", element_size)
        .staticmethod("element_size")
        .def("accessor", accessor)
        .def("nd", nd)
        .def("is_0_based", is_0_based)
        .def("origin", origin)
        .def("all", all)
        .def("last", last_0)
        .def("last", last_1)
        .def("is_padded", is_padded)
        .def("focus", focus_0)
        .def("focus", focus_1)
        .def("focus_size_1d", focus_size_1d)
        .def("is_trivial_1d", is_trivial_1d)
        .def("shift_origin", shift_origin)
        .def("id", id)
        .def("size", size)
        .def("__len__", size)
        .def("capacity", capacity)
        .def("__getitem__", getitem_1d, GetitemReturnValuePolicy())
        .def("__getitem__", getitem_1d_slice)
        .def("__getitem__", getitem_flex_grid, GetitemReturnValuePolicy())
        .def("__setitem__", setitem_1d)
        .def("__setitem__", setitem_flex_grid)
        .def("__delitem__", delitem_1d)
        .def("__delitem__", delitem_1d_slice)
        .def("front", front, GetitemReturnValuePolicy())
        .def("back", back, GetitemReturnValuePolicy())
        .def("fill", fill)
        .def("reserve", reserve)
        .def("deep_copy", deep_copy)
        .def("shallow_copy", shallow_copy)
        .def("as_1d", as_1d)
        .def("assign", assign)
        .def("pop_back", pop_back)
        .def("append", append)
        .def("insert", insert_i_x)
        .def("insert", insert_i_n_x)
        .def("resize", resize_1d_1)
        .def("resize", resize_1d_2)
        .def("resize", resize_flex_grid_1)
        .def("resize", resize_flex_grid_2)
        .def("reshape", reshape)
        .def("clear", clear)
        .def("extend", extend)
        .def("concatenate", concatenate)
        .def("reversed", reversed)
      ;
      select_wrappers<ElementType, f_t>::wrap(result);
      result
        .def("set_selected", set_selected_bool_a)
        .def("set_selected", set_selected_bool_s)
        .def("set_selected",
          (object(*)(
            object const&,
            af::const_ref<unsigned> const&,
            af::const_ref<e_t> const&)) set_selected_unsigned_a)
        .def("set_selected",
          (object(*)(
            object const&,
            af::const_ref<unsigned> const&,
            e_t const&)) set_selected_unsigned_s)
        .def("copy_selected",
          (object(*)(
            object const&,
            af::const_ref<unsigned> const&,
            af::const_ref<e_t> const&)) copy_selected_unsigned_a)
#if !defined(BOOST_ADAPTBX_TYPE_ID_SIZE_T_EQ_UNSIGNED)
        .def("set_selected",
          (object(*)(
            object const&,
            af::const_ref<std::size_t> const&,
            af::const_ref<e_t> const&)) set_selected_unsigned_a)
        .def("set_selected",
          (object(*)(
            object const&,
            af::const_ref<std::size_t> const&,
            e_t const&)) set_selected_unsigned_s)
        .def("copy_selected",
          (object(*)(
            object const&,
            af::const_ref<std::size_t> const&,
            af::const_ref<e_t> const&)) copy_selected_unsigned_a)
#endif
      ;
      return result;
    }

    static class_f_t
    ordered(std::string const& python_name,
            boost::python::object const& flex_root_scope)
    {
      {
        using namespace boost::python;
        scope local_scope(flex_root_scope);
        def("order", order_a_a);
        def("first_index", first_index_a_s);
        def("last_index", last_index_a_s);
      }
      return plain(python_name)
        .def("__eq__", eq_a_a)
        .def("__ne__", ne_a_a)
        .def("__eq__", eq_a_s)
        .def("__ne__", ne_a_s)
        .def("all_eq", all_eq_a_a)
        .def("all_ne", all_ne_a_a)
        .def("all_eq", all_eq_a_s)
        .def("all_ne", all_ne_a_s)
      ;
    }

    static class_f_t
    numeric_common(std::string const& python_name,
                   boost::python::object const& flex_root_scope)
    {
      {
        boost::python::scope local_scope(flex_root_scope);
        boost::python::def("sum", sum_a);
        boost::python::def("sum_sq", sum_sq_a);
        boost::python::def("product", product_a);
      }
      return plain(python_name)
        .def("count", count)
        .def("__neg__", neg_a)
        .def("__add__", add_a_a)
        .def("__sub__", sub_a_a)
        .def("__mul__", mul_a_a)
        .def("__div__", div_a_a)
        .def("__truediv__", div_a_a)
        .def("__add__", add_a_s)
        .def("__radd__", add_a_s)
        .def("__sub__", sub_a_s)
        .def("__rsub__", rsub_a_s)
        .def("__mul__", mul_a_s)
        .def("__rmul__", mul_a_s)
        .def("__div__", div_a_s)
        .def("__truediv__", div_a_s)
        .def("__rdiv__", rdiv_a_s)
        .def("__rtruediv__", rdiv_a_s)
        .def("__iadd__", iadd_a_a)
        .def("__isub__", isub_a_a)
        .def("__imul__", imul_a_a)
        .def("__idiv__", idiv_a_a)
        .def("__itruediv__", idiv_a_a)
        .def("__iadd__", iadd_a_s)
        .def("__isub__", isub_a_s)
        .def("__imul__", imul_a_s)
        .def("__idiv__", idiv_a_s)
        .def("__itruediv__", idiv_a_s)
        .def("__eq__", eq_a_a)
        .def("__ne__", ne_a_a)
        .def("__eq__", eq_a_s)
        .def("__ne__", ne_a_s)
        .def("all_eq", all_eq_a_a)
        .def("all_ne", all_ne_a_a)
        .def("all_eq", all_eq_a_s)
        .def("all_ne", all_ne_a_s)
      ;
    }

    static class_f_t
    numeric_no_pow(std::string const& python_name,
                   boost::python::object const& flex_root_scope)
    {
      {
        using namespace boost::python;
        scope local_scope(flex_root_scope);
        def("min_index", min_index_a);
        def("max_index", max_index_a);
        def("min", min_a);
        def("max", max_a);
        def("max_absolute", max_absolute_a);
        def("pow2", pow2_a);
        def("order", order_a_a);
        def("first_index", first_index_a_s);
        def("last_index", last_index_a_s);
      }
      return numeric_common(python_name, flex_root_scope)
        .def("as_double", as_double)
        .def("__lt__", lt_a_a)
        .def("__gt__", gt_a_a)
        .def("__le__", le_a_a)
        .def("__ge__", ge_a_a)
        .def("__lt__", lt_a_s)
        .def("__gt__", gt_a_s)
        .def("__le__", le_a_s)
        .def("__ge__", ge_a_s)
        .def("all_lt", all_lt_a_a)
        .def("all_gt", all_gt_a_a)
        .def("all_le", all_le_a_a)
        .def("all_ge", all_ge_a_a)
        .def("all_lt", all_lt_a_s)
        .def("all_gt", all_gt_a_s)
        .def("all_le", all_le_a_s)
        .def("all_ge", all_ge_a_s)
        .def("norm_inf", max_absolute_a)
      ;
    }

    static class_f_t
    numeric(std::string const& python_name,
            boost::python::object const& flex_root_scope)
    {
      {
        boost::python::scope local_scope(flex_root_scope);
        boost::python::def("abs", abs_a);
        boost::python::def("pow", pow_a_s);
        boost::python::def("fmod", fmod_a_s);
        boost::python::def("fmod_positive", fmod_positive_a_s);
        boost::python::def("atan2", atan2_a_a);
        boost::python::def("acos", acos_a);
        boost::python::def("cos", cos_a);
        boost::python::def("tan", tan_a);
        boost::python::def("asin", asin_a);
        boost::python::def("cosh", cosh_a);
        boost::python::def("tanh", tanh_a);
        boost::python::def("atan", atan_a);
        boost::python::def("exp", exp_a);
        boost::python::def("sin", sin_a);
        boost::python::def("fabs", fabs_a);
        boost::python::def("log", log_a);
        boost::python::def("sinh", sinh_a);
        boost::python::def("ceil", ceil_a);
        boost::python::def("floor", floor_a);
        boost::python::def("log10", log10_a);
        boost::python::def("sqrt", sqrt_a);
        boost::python::def("mean", mean_a);
        boost::python::def("mean_sq", mean_sq_a);
        boost::python::def("mean_weighted", mean_weighted_a_a);
        boost::python::def("mean_sq_weighted", mean_sq_weighted_a_a);
      }
      return numeric_no_pow(python_name, flex_root_scope)
        .def("norm", norm_a)
        .def("__pow__", pow_a_s)
      ;
    }

    static class_f_t
    integer(std::string const& python_name,
            boost::python::object const& flex_root_scope)
    {
      return numeric_no_pow(python_name, flex_root_scope)
        .def("__mod__", mod_a_a)
        .def("__mod__", mod_a_s)
        .def("__rmod__", rmod_a_s)
        .def("__imod__", imod_a_a)
        .def("__imod__", imod_a_s)
      ;
    }

    static class_f_t
    signed_integer(std::string const& python_name,
                   boost::python::object const& flex_root_scope)
    {
      {
        boost::python::scope local_scope(flex_root_scope);
        boost::python::def("abs", abs_a);
      }
      return integer(python_name, flex_root_scope);
    }
  };

}}} // namespace scitbx::af::boost_python

#endif // SCITBX_ARRAY_FAMILY_BOOST_PYTHON_FLEX_WRAPPER_H
