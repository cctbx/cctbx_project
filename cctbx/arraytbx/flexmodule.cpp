// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Aug: Created, based on sharedmodule.cpp, shared_bpl.h (rwgk)
 */

#include <boost/python/cross_module.hpp>
#include <cctbx/array_family/tiny_bpl.h>
#include <cctbx/array_family/small_bpl.h>
#include <cctbx/array_family/flex_types.h>
#include <cctbx/math/array_utils.h>

#include <cctbx/miller_bpl.h>
#include <cctbx/hendrickson_lattman_bpl.h>

#include <cctbx/sgtbx/matrix.h>
#include <cctbx/sftbx/xray_scatterer.h>

# include <cctbx/basic/from_bpl_import.h>

namespace cctbx { namespace af {

  boost::python::tuple flex_bool_getstate(
    flex_bool const& a);
  void flex_bool_setstate(
    flex_bool& a, boost::python::tuple state);
  boost::python::tuple flex_int_getstate(
    flex_int const& a);
  void flex_int_setstate(
    flex_int& a, boost::python::tuple state);
  boost::python::tuple flex_long_getstate(
    flex_long const& a);
  void flex_long_setstate(
    flex_long& a, boost::python::tuple state);
  boost::python::tuple flex_float_getstate(
    flex_float const& a);
  void flex_float_setstate(
    flex_float& a, boost::python::tuple state);
  boost::python::tuple flex_double_getstate(
    flex_double const& a);
  void flex_double_setstate(
    flex_double& a, boost::python::tuple state);
  boost::python::tuple flex_complex_double_getstate(
    flex_complex_double const& a);
  void flex_complex_double_setstate(
    flex_complex_double& a,
    boost::python::tuple state);
  boost::python::tuple flex_miller_index_getstate(
    versa<miller::Index, flex_grid<> > const& a);
  void flex_miller_index_setstate(
    versa<miller::Index, flex_grid<> >& a,
    boost::python::tuple state);
  boost::python::tuple flex_hendrickson_lattman_double_getstate(
    versa<hendrickson_lattman<double>, flex_grid<> > const& a);
  void flex_hendrickson_lattman_double_setstate(
    versa<hendrickson_lattman<double>, flex_grid<> >& a,
    boost::python::tuple state);
  boost::python::tuple flex_xray_scatterer_double_wk1995_getstate(
    versa<sftbx::XrayScatterer<double, eltbx::CAASF_WK1995>,
          flex_grid<> > const& a);
  void flex_xray_scatterer_double_wk1995_setstate(
    versa<sftbx::XrayScatterer<double, eltbx::CAASF_WK1995>,
          flex_grid<> >& a,
    boost::python::tuple state);

  template <typename ElementType>
  struct flex_pickle
  {
    template <typename ClassBuilderType>
    static void def(ClassBuilderType& class_bldr) {}
  };

  template <>
  struct flex_pickle<bool>
  {
    template <typename ClassBuilderType>
    static void def(ClassBuilderType& class_bldr)
    {
      class_bldr.def(flex_bool_getstate, "__getstate__");
      class_bldr.def(flex_bool_setstate, "__setstate__");
    }
  };

  template <>
  struct flex_pickle<int>
  {
    template <typename ClassBuilderType>
    static void def(ClassBuilderType& class_bldr)
    {
      class_bldr.def(flex_int_getstate, "__getstate__");
      class_bldr.def(flex_int_setstate, "__setstate__");
    }
  };

  template <>
  struct flex_pickle<long>
  {
    template <typename ClassBuilderType>
    static void def(ClassBuilderType& class_bldr)
    {
      class_bldr.def(flex_long_getstate, "__getstate__");
      class_bldr.def(flex_long_setstate, "__setstate__");
    }
  };

  template <>
  struct flex_pickle<float>
  {
    template <typename ClassBuilderType>
    static void def(ClassBuilderType& class_bldr)
    {
      class_bldr.def(flex_float_getstate, "__getstate__");
      class_bldr.def(flex_float_setstate, "__setstate__");
    }
  };

  template <>
  struct flex_pickle<double>
  {
    template <typename ClassBuilderType>
    static void def(ClassBuilderType& class_bldr)
    {
      class_bldr.def(flex_double_getstate, "__getstate__");
      class_bldr.def(flex_double_setstate, "__setstate__");
    }
  };

  template <>
  struct flex_pickle<std::complex<double> >
  {
    template <typename ClassBuilderType>
    static void def(ClassBuilderType& class_bldr)
    {
      class_bldr.def(flex_complex_double_getstate, "__getstate__");
      class_bldr.def(flex_complex_double_setstate, "__setstate__");
    }
  };

  template <>
  struct flex_pickle<miller::Index>
  {
    template <typename ClassBuilderType>
    static void def(ClassBuilderType& class_bldr)
    {
      class_bldr.def(flex_miller_index_getstate, "__getstate__");
      class_bldr.def(flex_miller_index_setstate, "__setstate__");
    }
  };

  template <>
  struct flex_pickle<hendrickson_lattman<double> >
  {
    template <typename ClassBuilderType>
    static void def(ClassBuilderType& class_bldr)
    {
      class_bldr.def(
        flex_hendrickson_lattman_double_getstate, "__getstate__");
      class_bldr.def(
        flex_hendrickson_lattman_double_setstate, "__setstate__");
    }
  };

  template <>
  struct flex_pickle<sftbx::XrayScatterer<double, eltbx::CAASF_WK1995> >
  {
    template <typename ClassBuilderType>
    static void def(ClassBuilderType& class_bldr)
    {
      class_bldr.def(
        flex_xray_scatterer_double_wk1995_getstate, "__getstate__");
      class_bldr.def(
        flex_xray_scatterer_double_wk1995_setstate, "__setstate__");
    }
  };

  struct flex_grid_wrappers
  {
    static
    flex_grid<>
    set_layout(
      flex_grid<>& fg,
      flex_grid_default_index_type const& layout)
    {
      return fg.set_layout(layout);
    }

    static
    flex_grid_default_index_type
    last_0(flex_grid<> const& fg)
    {
      return fg.last();
    }

    static
    flex_grid_default_index_type
    last_1(flex_grid<> const& fg, bool open_range)
    {
      return fg.last(open_range);
    }

    static
    tuple
    getinitargs(flex_grid<> const& fg)
    {
      bool open_range = true;
      tuple initargs(3);
      initargs.set_item(0, boost::python::make_ref(fg.origin()));
      initargs.set_item(1, boost::python::make_ref(fg.last(open_range)));
      initargs.set_item(2, boost::python::make_ref(open_range));
      return initargs;
    }

    static
    flex_grid_default_index_type
    getstate(flex_grid<> const& fg)
    {
      return fg.layout();
    }

    static
    void
    setstate(flex_grid<>& fg, flex_grid_default_index_type const& state)
    {
      fg.set_layout(state);
    }

  };

  template <typename ElementType>
  struct flex_items
  {
    flex_items() {}

    flex_items(versa<ElementType, flex_grid<> > const& data)
    : data_(data)
    {}

    std::size_t size() const { return data_.size(); }

    boost::python::tuple
    getitem(std::size_t i) const
    {
      if (i >= data_.size()) bpl_utils::raise_index_error();
      return boost::python::tuple(
        boost::python::make_ref(i), boost::python::make_ref(data_[i]));
    }

    versa<ElementType, flex_grid<> > data_;
  };

  template <typename ElementType>
  struct flex_wrapper : versa<ElementType, flex_grid<> >
  {
    typedef ElementType e_t;
    typedef versa<ElementType, flex_grid<> > f_t;
    typedef typename f_t::base_array_type base_array_type;

    typedef
    std::pair<
      boost::python::class_builder<f_t, flex_wrapper<ElementType> >,
      boost::python::class_builder<flex_items<ElementType> >
    >
    f_class_builders;

    // Tell the compiler how to convert a base class object to
    // this wrapper object.
    flex_wrapper(PyObject*, f_t const& v)
      : f_t(v)
    {}

    flex_wrapper(PyObject*)
      : f_t(flex_grid<>(0))
    {}

    flex_wrapper(PyObject*, flex_grid<> const& fg)
      : f_t(fg)
    {}

    flex_wrapper(PyObject*, flex_grid<> const& fg, ElementType const& x)
      : f_t(fg, x)
    {}

    flex_wrapper(PyObject*, std::size_t n)
      : f_t(flex_grid<>(n))
    {}

    flex_wrapper(PyObject*, std::size_t n, ElementType const& x)
      : f_t(flex_grid<>(n), x)
    {}

    flex_wrapper(PyObject*, boost::python::tuple tuple)
      : f_t(flex_grid<>(tuple.size()))
    {
      f_t::iterator a = this->begin();
      for(std::size_t i=0;i<tuple.size();i++)
        a[i] = BOOST_PYTHON_CONVERSION::from_python(
          tuple[i].get(), boost::python::type<e_t const&>());
    }

    static
    flex_grid<>
    accessor(f_t const& a) { return a.accessor(); }

    static
    std::size_t
    nd(f_t const& a) { return a.accessor().nd(); }

    static
    flex_grid_default_index_type
    origin(f_t const& a) { return a.accessor().origin(); }

    static
    flex_grid_default_index_type
    grid(f_t const& a) { return a.accessor().grid(); }

    static
    flex_grid_default_index_type
    last_0(f_t const& a) { return a.accessor().last(); }

    static
    flex_grid_default_index_type
    last_1(f_t const& a, bool open_range)
    {
      return a.accessor().last(open_range);
    }

    static
    flex_grid_default_index_type
    layout(f_t const& a) { return a.accessor().layout(); }

    static
    bool
    is_0_based(f_t const& a) { return a.accessor().is_0_based(); }

    static
    bool
    is_padded(f_t const& a) { return a.accessor().is_padded(); }

    static
    std::size_t
    id(f_t const& a) { return a.id(); }

    static
    std::size_t
    size(f_t const& a) { return a.size(); }

    static
    std::size_t
    capacity(f_t const& a) { return a.capacity(); }

    static
    e_t
    getitem_1d(f_t const& a, std::size_t i)
    {
      if (i >= a.size()) bpl_utils::raise_index_error();
      return a[i];
    }

    static
    e_t
    getitem_flex_grid(f_t const& a, flex_grid_default_index_type const& i)
    {
      if (!a.accessor().is_valid_index(i)) bpl_utils::raise_index_error();
      return a(i);
    }

    static
    void
    setitem_1d(f_t& a, std::size_t i, e_t const& x)
    {
      if (i >= a.size()) bpl_utils::raise_index_error();
      a[i] = x;
    }

    static
    void
    setitem_flex_grid(
      f_t& a, flex_grid_default_index_type const& i, e_t const& x)
    {
      if (!a.accessor().is_valid_index(i)) bpl_utils::raise_index_error();
      a(i) = x;
    }

    static
    e_t
    front(f_t const& a)
    {
      if (a.size() == 0) bpl_utils::raise_index_error();
      return a.front();
    }

    static
    e_t
    back(f_t const& a)
    {
      if (a.size() == 0) bpl_utils::raise_index_error();
      return a.back();
    }

    static
    void
    fill(f_t& a, e_t const& x)
    {
      a.fill(x);
    }

    static
    void
    reserve(f_t& a, std::size_t sz)
    {
      a.reserve(sz);
    }

    static
    f_t
    deep_copy(f_t const& a)
    {
      return a.deep_copy();
    }

    static
    f_t
    shallow_copy(f_t const& a)
    {
      return a;
    }

    static
    f_t
    as_1d(f_t const& a)
    {
      flex_grid_default_index_type grid;
      grid.push_back(a.size());
      return f_t(a, flex_grid<>(grid));
    }

    static
    void
    assign(f_t& a, std::size_t sz, e_t const& x)
    {
      base_array_type b = bpl_utils::as_base_array(a);
      b.assign(sz, x);
      a.resize(flex_grid<>(b.size()));
    }

    static
    void push_back(f_t& a, e_t const& x)
    {
      base_array_type b = bpl_utils::as_base_array(a);
      b.push_back(x);
      a.resize(flex_grid<>(b.size()));
    }

    static
    void
    pop_back(f_t& a)
    {
      base_array_type b = bpl_utils::as_base_array(a);
      if (b.size() == 0) bpl_utils::raise_index_error();
      b.pop_back();
      a.resize(flex_grid<>(b.size()));
    }

    static
    void
    insert_i_x(f_t& a, std::size_t i, e_t const& x)
    {
      base_array_type b = bpl_utils::as_base_array(a);
      if (i >= b.size()) bpl_utils::raise_index_error();
      b.insert(&b[i], x);
      a.resize(flex_grid<>(b.size()));
    }

    static
    void
    insert_i_n_x(f_t& a, std::size_t i, std::size_t n, e_t const& x)
    {
      base_array_type b = bpl_utils::as_base_array(a);
      if (i >= b.size()) bpl_utils::raise_index_error();
      b.insert(&b[i], n, x);
      a.resize(flex_grid<>(b.size()));
    }

    static
    void
    erase_i(f_t& a, std::size_t i)
    {
      base_array_type b = bpl_utils::as_base_array(a);
      if (i >= b.size()) bpl_utils::raise_index_error();
      b.erase(&b[i]);
      a.resize(flex_grid<>(b.size()));
    }

    static
    void erase_i_j(f_t& a, std::size_t i, std::size_t j)
    {
      base_array_type b = bpl_utils::as_base_array(a);
      if (i >= b.size()) bpl_utils::raise_index_error();
      if (j >= b.size()) bpl_utils::raise_index_error();
      b.erase(&b[i], &b[j]);
      a.resize(flex_grid<>(b.size()));
    }

    static
    void
    resize_1d_1(f_t& a, std::size_t sz)
    {
      base_array_type b = bpl_utils::as_base_array(a);
      b.resize(sz);
      a.resize(flex_grid<>(b.size()));
    }

    static
    void
    resize_1d_2(f_t& a, std::size_t sz, e_t const& x)
    {
      base_array_type b = bpl_utils::as_base_array(a);
      b.resize(sz, x);
      a.resize(flex_grid<>(b.size()));
    }

    static
    void
    resize_flex_grid_1(f_t& a, flex_grid<> const& grid)
    {
      a.resize(grid);
    }

    static
    void
    resize_flex_grid_2(f_t& a, flex_grid<> const& grid, e_t const& x)
    {
      a.resize(grid, x);
    }

    static
    void
    clear(f_t& a)
    {
      base_array_type b = bpl_utils::as_base_array(a);
      b.clear();
      a.resize(flex_grid<>(b.size()));
    }

    static
    void
    append(f_t& a, f_t const& other)
    {
      base_array_type b = bpl_utils::as_base_array(a);
      bpl_utils::assert_0_based_1d(other.accessor());
      b.insert(b.end(), other.begin(), other.end());
      a.resize(flex_grid<>(b.size()));
    }

    static
    boost::python::ref
    indices(f_t const& a)
    {
      return boost::python::ref(PyRange_New(0, a.size(), 1, 1));
    }

    static
    flex_items<e_t>
    items(f_t const& a)
    {
      return flex_items<e_t>(a);
    }

    static
    f_t
    select(f_t const& a, flex_bool const& flags)
    {
      bpl_utils::assert_0_based_1d(a.accessor());
      bpl_utils::assert_0_based_1d(flags.accessor());
      if (a.size() != flags.size()) bpl_utils::raise_incompatible_arrays();
      std::size_t n = 0;
      std::size_t i;
      for(i=0;i<flags.size();i++) if (flags[i]) n++;
      base_array_type result;
      result.reserve(n);
      for(i=0;i<flags.size();i++) if (flags[i]) result.push_back(a[i]);
      return f_t(result, flex_grid<>(result.size()));
    }

    static
    f_t
    shuffle(f_t const& a, flex_size_t const& permutation)
    {
      bpl_utils::assert_0_based_1d(a.accessor());
      bpl_utils::assert_0_based_1d(permutation.accessor());
      if (a.size() != permutation.size()) {
        bpl_utils::raise_incompatible_arrays();
      }
      base_array_type result;
      if (a.size()) {
        result.resize(a.size(), a[0]); // avoid requirement that e_t is
        for(std::size_t i=1;i<a.size();i++) {  // default constructible
          std::size_t j = permutation[i];
          cctbx_assert(j < a.size());
          result[j] = a[i];
        }
      }
      return f_t(result, a.accessor());
    }

    static
    flex_bool
    invert_a(flex_bool const& a)
    {
      shared_plain<bool> result;
      result.reserve(a.size());
      for(std::size_t i=0;i<a.size();i++) result.push_back(!a[i]);
      return flex_bool(result, a.accessor());
    }

    static
    flex_bool
    and_a_a(flex_bool const& a1, flex_bool const& a2)
    {
      if (a1.accessor() != a2.accessor()) {
        bpl_utils::raise_incompatible_arrays();
      }
      shared_plain<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.push_back(a1[i] && a2[i]);
      return flex_bool(result, a1.accessor());
    }

    static
    flex_bool
    or_a_a(flex_bool const& a1, flex_bool const& a2)
    {
      if (a1.accessor() != a2.accessor()) {
        bpl_utils::raise_incompatible_arrays();
      }
      shared_plain<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.push_back(a1[i] || a2[i]);
      return flex_bool(result, a1.accessor());
    }

    static
    flex_bool
    iand_a_a(flex_bool a1, flex_bool const& a2)
    {
      if (a1.accessor() != a2.accessor()) {
        bpl_utils::raise_incompatible_arrays();
      }
      for(std::size_t i=0;i<a1.size();i++) if(!a2[i]) a1[i] = false;
      return a1;
    }

    static
    flex_bool
    ior_a_a(flex_bool a1, flex_bool const& a2)
    {
      if (a1.accessor() != a2.accessor()) {
        bpl_utils::raise_incompatible_arrays();
      }
      for(std::size_t i=0;i<a1.size();i++) if(a2[i]) a1[i] = true;
      return a1;
    }

    static
    flex_bool
    iand_a_s(flex_bool a1, bool a2)
    {
      if (!a2) std::fill(a1.begin(), a1.end(), false);
      return a1;
    }

    static
    flex_bool
    ior_a_s(flex_bool a1, bool a2)
    {
      if (a2) std::fill(a1.begin(), a1.end(), true);
      return a1;
    }

    static
    std::size_t
    count(flex_bool const& a1, bool a2)
    {
      std::size_t result = 0;
      for(std::size_t i=0;i<a1.size();i++) if(a1[i] == a2) result++;
      return result;
    }

    static
    flex_double
    as_double(f_t const& a)
    {
      shared_plain<double> result;
      result.reserve(a.size());
      for(std::size_t i=0;i<a.size();i++) result.push_back(a[i]);
      return flex_double(result, a.accessor());
    }

    static
    f_t
    neg_a(f_t const& a)
    {
      base_array_type result;
      result.reserve(a.size());
      for(std::size_t i=0;i<a.size();i++) result.push_back(-a[i]);
      return f_t(result, a.accessor());
    }

    static
    f_t
    add_a_a(f_t const& a1, f_t const& a2)
    {
      if (a1.accessor() != a2.accessor()) {
        bpl_utils::raise_incompatible_arrays();
      }
      base_array_type result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.push_back(a1[i] + a2[i]);
      return f_t(result, a1.accessor());
    }

    static
    f_t
    sub_a_a(f_t const& a1, f_t const& a2)
    {
      if (a1.accessor() != a2.accessor()) {
        bpl_utils::raise_incompatible_arrays();
      }
      base_array_type result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.push_back(a1[i] - a2[i]);
      return f_t(result, a1.accessor());
    }

    static
    f_t
    mul_a_a(f_t const& a1, f_t const& a2)
    {
      if (a1.accessor() != a2.accessor()) {
        bpl_utils::raise_incompatible_arrays();
      }
      base_array_type result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.push_back(a1[i] * a2[i]);
      return f_t(result, a1.accessor());
    }

    static
    f_t
    div_a_a(f_t const& a1, f_t const& a2)
    {
      if (a1.accessor() != a2.accessor()) {
        bpl_utils::raise_incompatible_arrays();
      }
      base_array_type result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.push_back(a1[i] / a2[i]);
      return f_t(result, a1.accessor());
    }

    static
    f_t
    mod_a_a(f_t const& a1, f_t const& a2)
    {
      if (a1.accessor() != a2.accessor()) {
        bpl_utils::raise_incompatible_arrays();
      }
      base_array_type result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.push_back(a1[i] % a2[i]);
      return f_t(result, a1.accessor());
    }

    static
    f_t
    add_a_s(f_t const& a1, e_t const& a2)
    {
      base_array_type result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.push_back(a1[i] + a2);
      return f_t(result, a1.accessor());
    }

    static
    f_t
    sub_a_s(f_t const& a1, e_t const& a2)
    {
      base_array_type result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.push_back(a1[i] - a2);
      return f_t(result, a1.accessor());
    }

    static
    f_t
    mul_a_s(f_t const& a1, e_t const& a2)
    {
      base_array_type result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.push_back(a1[i] * a2);
      return f_t(result, a1.accessor());
    }

    static
    f_t
    div_a_s(f_t const& a1, e_t const& a2)
    {
      base_array_type result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.push_back(a1[i] / a2);
      return f_t(result, a1.accessor());
    }

    static
    f_t
    mod_a_s(f_t const& a1, e_t const& a2)
    {
      base_array_type result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.push_back(a1[i] % a2);
      return f_t(result, a1.accessor());
    }

    static
    f_t
    iadd_a_s(f_t& a1, e_t const& a2)
    {
      for(std::size_t i=0;i<a1.size();i++) a1[i] += a2;
      return a1;
    }

    static
    f_t
    isub_a_s(f_t& a1, e_t const& a2)
    {
      for(std::size_t i=0;i<a1.size();i++) a1[i] -= a2;
      return a1;
    }

    static
    f_t
    imul_a_s(f_t& a1, e_t const& a2)
    {
      for(std::size_t i=0;i<a1.size();i++) a1[i] *= a2;
      return a1;
    }

    static
    f_t
    idiv_a_s(f_t& a1, e_t const& a2)
    {
      for(std::size_t i=0;i<a1.size();i++) a1[i] /= a2;
      return a1;
    }

    static
    f_t
    imod_a_s(f_t& a1, e_t const& a2)
    {
      for(std::size_t i=0;i<a1.size();i++) a1[i] %= a2;
      return a1;
    }

    static
    int
    cmp_a_a(f_t const& a1, f_t const& a2)
    {
      return af::cmp(a1.const_ref(), a2.const_ref());
    }

    static
    int
    cmp_a_s(f_t const& a1, e_t const& a2)
    {
      return af::cmp(a1.const_ref(), a2);
    }

    static
    flex_bool
    eq_a_a(f_t const& a1, f_t const& a2)
    {
      if (a1.accessor() != a2.accessor()) {
        bpl_utils::raise_incompatible_arrays();
      }
      shared_plain<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.push_back(a1[i] == a2[i]);
      return flex_bool(result, a1.accessor());
    }

    static
    flex_bool
    ne_a_a(f_t const& a1, f_t const& a2)
    {
      if (a1.accessor() != a2.accessor()) {
        bpl_utils::raise_incompatible_arrays();
      }
      shared_plain<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.push_back(a1[i] != a2[i]);
      return flex_bool(result, a1.accessor());
    }

    static
    flex_bool
    lt_a_a(f_t const& a1, f_t const& a2)
    {
      if (a1.accessor() != a2.accessor()) {
        bpl_utils::raise_incompatible_arrays();
      }
      shared_plain<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.push_back(a1[i] < a2[i]);
      return flex_bool(result, a1.accessor());
    }

    static
    flex_bool
    gt_a_a(f_t const& a1, f_t const& a2)
    {
      if (a1.accessor() != a2.accessor()) {
        bpl_utils::raise_incompatible_arrays();
      }
      shared_plain<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.push_back(a1[i] > a2[i]);
      return flex_bool(result, a1.accessor());
    }

    static
    flex_bool
    le_a_a(f_t const& a1, f_t const& a2)
    {
      if (a1.accessor() != a2.accessor()) {
        bpl_utils::raise_incompatible_arrays();
      }
      shared_plain<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.push_back(a1[i] <= a2[i]);
      return flex_bool(result, a1.accessor());
    }

    static
    flex_bool
    ge_a_a(f_t const& a1, f_t const& a2)
    {
      if (a1.accessor() != a2.accessor()) {
        bpl_utils::raise_incompatible_arrays();
      }
      shared_plain<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.push_back(a1[i] >= a2[i]);
      return flex_bool(result, a1.accessor());
    }

    static
    flex_bool
    eq_a_s(f_t const& a1, e_t const& a2)
    {
      shared_plain<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.push_back(a1[i] == a2);
      return flex_bool(result, a1.accessor());
    }

    static
    flex_bool
    ne_a_s(f_t const& a1, e_t const& a2)
    {
      shared_plain<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.push_back(a1[i] != a2);
      return flex_bool(result, a1.accessor());
    }

    static
    flex_bool
    lt_a_s(f_t const& a1, e_t const& a2)
    {
      shared_plain<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.push_back(a1[i] < a2);
      return flex_bool(result, a1.accessor());
    }

    static
    flex_bool
    gt_a_s(f_t const& a1, e_t const& a2)
    {
      shared_plain<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.push_back(a1[i] > a2);
      return flex_bool(result, a1.accessor());
    }

    static
    flex_bool
    le_a_s(f_t const& a1, e_t const& a2)
    {
      shared_plain<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.push_back(a1[i] <= a2);
      return flex_bool(result, a1.accessor());
    }

    static
    flex_bool
    ge_a_s(f_t const& a1, e_t const& a2)
    {
      shared_plain<bool> result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) result.push_back(a1[i] >= a2);
      return flex_bool(result, a1.accessor());
    }

    static
    f_t
    abs_a(f_t const& a)
    {
      base_array_type result;
      result.reserve(a.size());
      for(std::size_t i=0;i<a.size();i++) result.append(math::abs(a[i]));
      return f_t(result, a.accessor());
    }

    static
    f_t
    fmod_a_s(f_t const& a1, e_t const& a2)
    {
      base_array_type result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) {
        result.append(std::fmod(a1[i], a2));
      }
      return f_t(result, a1.accessor());
    }

    static
    f_t
    pow_a_s(f_t const& a, e_t const& exponent)
    {
      base_array_type result;
      result.reserve(a.size());
      for(std::size_t i=0;i<a.size();i++) {
        result.append(std::pow(a[i], exponent));
      }
      return f_t(result, a.accessor());
    }

    static
    f_t
    atan2_a_a(f_t const& a1, f_t const& a2)
    {
      if (a1.accessor() != a2.accessor()) {
        bpl_utils::raise_incompatible_arrays();
      }
      base_array_type result;
      result.reserve(a1.size());
      for(std::size_t i=0;i<a1.size();i++) {
        result.append(std::atan2(a1[i], a2[i]));
      }
      return f_t(result, a1.accessor());
    }

#define CCTBX_ARRAY_FAMILY_SHARED_BPL_CMATH_1ARG(func) \
    static \
    f_t \
    func ## _a(f_t const& a) \
    { \
      base_array_type result; \
      result.reserve(a.size()); \
      for(std::size_t i=0;i<a.size();i++) { \
        result.append(std::func(a[i])); \
      } \
      return f_t(result, a.accessor()); \
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
    max_index_a(f_t const& a) { return af::max_index(a.const_ref()); }
    static std::size_t
    min_index_a(f_t const& a) { return af::min_index(a.const_ref()); }
    static e_t max_a(f_t const& a) { return af::max(a.const_ref()); }
    static e_t min_a(f_t const& a) { return af::min(a.const_ref()); }
    static e_t sum_a(f_t const& a) { return af::sum(a.const_ref()); }
    static e_t sum_sq_a(f_t const& a) { return af::sum_sq(a.const_ref()); }
    static e_t product_a(f_t const& a) { return af::product(a.const_ref()); }
    static e_t mean_a(f_t const& a) { return af::mean(a.const_ref()); }
    static e_t mean_sq_a(f_t const& a) { return af::mean_sq(a.const_ref()); }

    static
    e_t
    mean_weighted_a_a(f_t const& a1, f_t const& a2)
    {
      return af::mean_weighted(a1.const_ref(), a2.const_ref());
    }

    static
    e_t
    mean_sq_weighted_a_a(f_t const& a1, f_t const& a2)
    {
      return af::mean_sq_weighted(a1.const_ref(), a2.const_ref());
    }

    static
    versa<double, flex_grid<> >
    real_complex(versa<std::complex<double>, flex_grid<> > const& a)
    {
      shared_plain<double> result;
      result.reserve(a.size());
      for(std::size_t i=0;i<a.size();i++) {
        result.push_back(std::real(a[i]));
      }
      return versa<double, flex_grid<> >(result, a.accessor());
    }

    static
    versa<double, flex_grid<> >
    imag_complex(versa<std::complex<double>, flex_grid<> > const& a)
    {
      shared_plain<double> result;
      result.reserve(a.size());
      for(std::size_t i=0;i<a.size();i++) {
        result.push_back(std::imag(a[i]));
      }
      return versa<double, flex_grid<> >(result, a.accessor());
    }

    static
    versa<double, flex_grid<> >
    abs_complex(versa<std::complex<double>, flex_grid<> > const& a)
    {
      shared_plain<double> result;
      result.reserve(a.size());
      for(std::size_t i=0;i<a.size();i++) {
        result.push_back(std::abs(a[i]));
      }
      return versa<double, flex_grid<> >(result, a.accessor());
    }

    static
    versa<double, flex_grid<> >
    arg_complex_2(versa<std::complex<double>, flex_grid<> > const& a, bool deg)
    {
      shared_plain<double> result;
      result.reserve(a.size());
      for(std::size_t i=0;i<a.size();i++) {
        result.push_back(std::arg(a[i]));
        if (deg) result[i] /= constants::pi_180;
      }
      return versa<double, flex_grid<> >(result, a.accessor());
    }

    static
    versa<double, flex_grid<> >
    arg_complex_1(versa<std::complex<double>, flex_grid<> > const& a)
    {
      return arg_complex_2(a, false);
    }

    static
    versa<double, flex_grid<> >
    norm_complex(versa<std::complex<double>, flex_grid<> > const& a)
    {
      shared_plain<double> result;
      result.reserve(a.size());
      for(std::size_t i=0;i<a.size();i++) {
        result.push_back(std::norm(a[i]));
      }
      return versa<double, flex_grid<> >(result, a.accessor());
    }

    static
    versa<std::complex<double>, flex_grid<> >
    polar_complex_3(
      versa<double, flex_grid<> > const& rho,
      versa<double, flex_grid<> > const& theta,
      bool deg)
    {
      if (rho.accessor() != theta.accessor()) {
        bpl_utils::raise_incompatible_arrays();
      }
      shared_plain<std::complex<double> > result;
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
      return versa<std::complex<double>, flex_grid<> >(result, rho.accessor());
    }

    static
    versa<std::complex<double>, flex_grid<> >
    polar_complex_2(
      versa<double, flex_grid<> > const& rho,
      versa<double, flex_grid<> > const& theta)
    {
      return polar_complex_3(rho, theta, false);
    }

    static
    f_class_builders
    plain(
      boost::python::module_builder& bpl_module,
      std::string const& python_name)
    {
      class_builder<flex_items<ElementType> >
      py_flex_items(bpl_module, (python_name+"_items").c_str());

      class_builder<f_t, flex_wrapper<ElementType> >
      py_flex(bpl_module, python_name.c_str());
      boost::python::export_converters(py_flex);

      py_flex_items.def(constructor<>());
      py_flex_items.def(constructor<f_t const&>());
      py_flex_items.def(&flex_items<ElementType>::size, "__len__");
      py_flex_items.def(&flex_items<ElementType>::getitem, "__getitem__");

      py_flex.def(constructor<>());
      py_flex.def(constructor<flex_grid<> const&>());
      py_flex.def(constructor<flex_grid<> const&, ElementType const&>());
      py_flex.def(constructor<std::size_t>());
      py_flex.def(constructor<std::size_t, ElementType const&>());
      py_flex.def(constructor<tuple>());
      py_flex.def(accessor, "accessor");
      py_flex.def(nd, "nd");
      py_flex.def(origin, "origin");
      py_flex.def(grid, "grid");
      py_flex.def(last_0, "last");
      py_flex.def(last_1, "last");
      py_flex.def(layout, "layout");
      py_flex.def(is_0_based, "is_0_based");
      py_flex.def(is_padded, "is_padded");
      py_flex.def(id, "id");
      py_flex.def(size, "size");
      py_flex.def(size, "__len__");
      py_flex.def(capacity, "capacity");
      py_flex.def(getitem_1d, "__getitem__");
      py_flex.def(getitem_flex_grid, "__getitem__");
      py_flex.def(setitem_1d, "__setitem__");
      py_flex.def(setitem_flex_grid, "__setitem__");
      py_flex.def(front, "front");
      py_flex.def(back, "back");
      py_flex.def(fill, "fill");
      py_flex.def(deep_copy, "deep_copy");
      py_flex.def(shallow_copy, "shallow_copy");
      py_flex.def(as_1d, "as_1d");
      py_flex.def(assign, "assign");
      py_flex.def(push_back, "push_back");
      py_flex.def(push_back, "append");
      py_flex.def(pop_back, "pop_back");
      py_flex.def(insert_i_x, "insert");
      py_flex.def(insert_i_n_x, "insert");
      py_flex.def(erase_i, "erase");
      py_flex.def(erase_i_j, "erase");
      py_flex.def(resize_1d_1, "resize");
      py_flex.def(resize_1d_2, "resize");
      py_flex.def(resize_flex_grid_1, "resize");
      py_flex.def(resize_flex_grid_2, "resize");
      py_flex.def(clear, "clear");
      py_flex.def(append, "append");
      py_flex.def(indices, "indices");
      py_flex.def(items, "items");
      py_flex.def(select, "select");
      py_flex.def(shuffle, "shuffle");

      flex_pickle<e_t>::def(py_flex);

      return std::make_pair(py_flex, py_flex_items);
    }

    static
    f_class_builders
    cmp_comparable(
      boost::python::module_builder& bpl_module,
      std::string const& python_name)
    {
      f_class_builders class_blds = plain(bpl_module, python_name);
      class_blds.first.def(cmp_a_a, "__cmp__");
      class_blds.first.def(cmp_a_s, "__cmp__");
      return class_blds;
    }

    static
    f_class_builders
    logical(
      boost::python::module_builder& bpl_module,
      std::string const& python_name)
    {
      f_class_builders class_blds = cmp_comparable(bpl_module, python_name);
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
    f_class_builders
    numeric_common(
      boost::python::module_builder& bpl_module,
      std::string const& python_name)
    {
      f_class_builders class_blds = plain(bpl_module, python_name);
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
    f_class_builders
    numeric_no_pow(
      boost::python::module_builder& bpl_module,
      std::string const& python_name)
    {
      f_class_builders class_blds = numeric_common(bpl_module, python_name);
      class_blds.first.def(as_double, "as_double");
      class_blds.first.def(cmp_a_a, "__cmp__");
      class_blds.first.def(cmp_a_s, "cmp"); // XXX __cmp__ did not work
      class_blds.first.def(eq_a_a, "__eq__");
      class_blds.first.def(ne_a_a, "__ne__");
      class_blds.first.def(eq_a_s, "__eq__");
      class_blds.first.def(ne_a_s, "__ne__");
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
    f_class_builders
    numeric(
      boost::python::module_builder& bpl_module,
      std::string const& python_name)
    {
      f_class_builders class_blds = numeric_no_pow(bpl_module, python_name);
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
    f_class_builders
    integer(
      boost::python::module_builder& bpl_module,
      std::string const& python_name)
    {
      f_class_builders class_blds = numeric_no_pow(bpl_module, python_name);
      class_blds.first.def(mod_a_a, "__mod__");
      class_blds.first.def(mod_a_s, "mod");
      class_blds.first.def(imod_a_s, "__imod__");
      return class_blds;
    }

    static
    f_class_builders
    complex(
      boost::python::module_builder& bpl_module,
      std::string const& python_name)
    {
      f_class_builders class_blds = numeric_common(bpl_module, python_name);
      class_blds.first.def(eq_a_a, "__eq__");
      class_blds.first.def(ne_a_a, "__ne__");
      class_blds.first.def(eq_a_s, "__eq__");
      class_blds.first.def(ne_a_s, "__ne__");
      bpl_module.def(real_complex, "real");
      bpl_module.def(imag_complex, "imag");
      bpl_module.def(abs_complex, "abs");
      bpl_module.def(arg_complex_2, "arg");
      bpl_module.def(arg_complex_1, "arg");
      bpl_module.def(norm_complex, "norm");
      bpl_module.def(polar_complex_3, "polar");
      bpl_module.def(polar_complex_2, "polar");
      return class_blds;
    }
  };

  void init_module(python::module_builder& this_module)
  {
    const std::string Revision = "$Revision$";
    this_module.add(python::ref(to_python(
        Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

    python::import_converters<cctbx::sgtbx::RTMx>
    py_RTMx("cctbx_boost.sgtbx", "RTMx");

    typedef cctbx::sftbx::XrayScatterer<
      double, cctbx::eltbx::CAASF_WK1995> XrayScatterer;
    python::import_converters<XrayScatterer>
    py_XrayScatterer("cctbx_boost.sftbx", "XrayScatterer");

    class_builder<flex_grid<> > py_flex_grid(this_module, "grid");

    py_flex_grid.def(constructor<>());
    py_flex_grid.def(constructor<
      flex_grid_default_index_type const&>());
    py_flex_grid.def(constructor<
      flex_grid_default_index_type const&,
      flex_grid_default_index_type const&>());
    py_flex_grid.def(constructor<
      flex_grid_default_index_type const&,
      flex_grid_default_index_type const&,
      bool>());
    py_flex_grid.def(flex_grid_wrappers::set_layout, "set_layout");
    py_flex_grid.def(&flex_grid<>::nd, "nd");
    py_flex_grid.def(&flex_grid<>::size1d, "size1d");
    py_flex_grid.def(&flex_grid<>::origin, "origin");
    py_flex_grid.def(&flex_grid<>::grid, "grid");
    py_flex_grid.def(flex_grid_wrappers::last_0, "last");
    py_flex_grid.def(flex_grid_wrappers::last_1, "last");
    py_flex_grid.def(&flex_grid<>::layout, "layout");
    py_flex_grid.def(&flex_grid<>::is_0_based, "is_0_based");
    py_flex_grid.def(&flex_grid<>::is_padded, "is_padded");
    py_flex_grid.def(&flex_grid<>::operator(), "__call__");
    py_flex_grid.def(&flex_grid<>::is_valid_index, "is_valid_index");
    py_flex_grid.def(&flex_grid<>::operator==, "__eq__");
    py_flex_grid.def(&flex_grid<>::operator!=, "__ne__");
    py_flex_grid.def(flex_grid_wrappers::getinitargs, "__getinitargs__");
    py_flex_grid.def(flex_grid_wrappers::getstate, "__getstate__");
    py_flex_grid.def(flex_grid_wrappers::setstate, "__setstate__");

#define WRAP_PLAIN(python_name, element_type) \
    flex_wrapper<element_type >::plain(this_module, python_name)
#define WRAP_CMP_COMPARABLE(python_name, element_type) \
    flex_wrapper<element_type >::cmp_comparable(this_module, python_name)
#define WRAP_INTEGER(python_name, element_type) \
    flex_wrapper<element_type >::integer(this_module, python_name)
#define WRAP_NUMERIC(python_name, element_type) \
    flex_wrapper<element_type >::numeric(this_module, python_name)
#define WRAP_COMPLEX(python_name, element_type) \
    flex_wrapper<element_type >::complex(this_module, python_name)

    // bool is wrapped here to enable boolean operators for numeric types
    flex_wrapper<bool>::logical(this_module, "bool");

    // size_t is wrapped here to enable .shuffle() for the other types
    typedef std::size_t size_t;
    WRAP_CMP_COMPARABLE("size_t", size_t);

    // double is wrapped here to enable .as_double() for the other types
    WRAP_NUMERIC("double", double);

#ifndef FAST_COMPILE
    WRAP_INTEGER("int", int);
    WRAP_NUMERIC("float", float);
    WRAP_COMPLEX("complex_double", std::complex<double>);

    WRAP_INTEGER("long", long);
    WRAP_CMP_COMPARABLE("std_string", std::string);

    WRAP_CMP_COMPARABLE("miller_Index", cctbx::miller::Index);
    WRAP_PLAIN("hendrickson_lattman", cctbx::hendrickson_lattman<double>);
    WRAP_PLAIN("RTMx", cctbx::sgtbx::RTMx);
    WRAP_PLAIN("XrayScatterer", XrayScatterer);
    WRAP_PLAIN("double3", cctbx::af::double3);

    typedef cctbx::af::tiny<size_t, 2> tiny_size_t_2;
    WRAP_PLAIN("tiny_size_t_2", tiny_size_t_2);
#endif
  }

  // hook for flex_picklers.cpp
  boost::python::ref
  make_ref_flex_grid(flex_grid<> const& fg)
  {
    return boost::python::make_ref(fg);
  }

  // hook for flex_picklers.cpp
  flex_grid<>
  flex_grid_from_python(PyObject* obj)
  {
    return BOOST_PYTHON_CONVERSION::from_python(
      obj, boost::python::type<flex_grid<> const&>());
  }

}} // namespace cctbx::af

BOOST_PYTHON_MODULE_INIT(flex)
{
  boost::python::module_builder this_module("flex");
  cctbx::af::init_module(this_module);
}
