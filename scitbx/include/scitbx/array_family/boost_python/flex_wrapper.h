/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copied from cctbx/arraytbx/flexmodule.cpp (rwgk)
     2002 Aug: Created, based on sharedmodule.cpp, shared_bpl.h (rwgk)
 */

#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_FLEX_WRAPPER_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_FLEX_WRAPPER_H

#include <scitbx/constants.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/versa_reductions.h>
#include <scitbx/array_family/versa_algebra.h>
#include <scitbx/boost_python/utils.h>
#include <scitbx/boost_python/container_conversions.h>
#include <scitbx/array_family/boost_python/shared_flex_conversions.h>
#include <scitbx/array_family/boost_python/ref_flex_conversions.h>
#include <scitbx/array_family/boost_python/utils.h>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_non_const_reference.hpp>

namespace scitbx { namespace af { namespace boost_python {

  using scitbx::boost_python::positive_getitem_index;

  template <typename ElementType>
  struct flex_items
  {
    flex_items() {}

    flex_items(versa<ElementType, flex_grid<> > const& data)
    : data_(data)
    {}

    std::size_t size() const { return data_.size(); }

    boost::python::tuple
    getitem(long i) const
    {
      std::size_t j = positive_getitem_index(i, data_.size());
      return boost::python::make_tuple(j, data_[j]);
    }

    versa<ElementType, flex_grid<> > data_;
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

    flex_wrapper(PyObject*, shared_plain<ElementType> const& a)
      : f_t(a, flex_grid<>(a.size()))
    {}

    static flex_grid<>
    accessor(f_t const& a) { return a.accessor(); }

    static std::size_t
    nd(f_t const& a) { return a.accessor().nd(); }

    static flex_grid_default_index_type
    origin(f_t const& a) { return a.accessor().origin(); }

    static flex_grid_default_index_type
    grid(f_t const& a) { return a.accessor().grid(); }

    static flex_grid_default_index_type
    last_0(f_t const& a) { return a.accessor().last(); }

    static flex_grid_default_index_type
    last_1(f_t const& a, bool open_range)
    {
      return a.accessor().last(open_range);
    }

    static flex_grid_default_index_type
    layout(f_t const& a) { return a.accessor().layout(); }

    static bool
    is_0_based(f_t const& a) { return a.accessor().is_0_based(); }

    static bool
    is_padded(f_t const& a) { return a.accessor().is_padded(); }

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
      return f_t(a, flex_grid<>(a.size()));
    }

    static void
    assign(f_t& a, std::size_t sz, e_t const& x)
    {
      base_array_type b = flex_as_base_array(a);
      b.assign(sz, x);
      a.resize(flex_grid<>(b.size()));
    }

    static void
    push_back(f_t& a, e_t const& x)
    {
      base_array_type b = flex_as_base_array(a);
      b.push_back(x);
      a.resize(flex_grid<>(b.size()));
    }

    static void
    pop_back(f_t& a)
    {
      base_array_type b = flex_as_base_array(a);
      if (b.size() == 0) scitbx::boost_python::raise_index_error();
      b.pop_back();
      a.resize(flex_grid<>(b.size()));
    }

    static void
    insert_i_x(f_t& a, std::size_t i, e_t const& x)
    {
      base_array_type b = flex_as_base_array(a);
      if (i >= b.size()) scitbx::boost_python::raise_index_error();
      b.insert(&b[i], x);
      a.resize(flex_grid<>(b.size()));
    }

    static void
    insert_i_n_x(f_t& a, std::size_t i, std::size_t n, e_t const& x)
    {
      base_array_type b = flex_as_base_array(a);
      if (i >= b.size()) scitbx::boost_python::raise_index_error();
      b.insert(&b[i], n, x);
      a.resize(flex_grid<>(b.size()));
    }

    static void
    erase_i(f_t& a, std::size_t i)
    {
      base_array_type b = flex_as_base_array(a);
      if (i >= b.size()) scitbx::boost_python::raise_index_error();
      b.erase(&b[i]);
      a.resize(flex_grid<>(b.size()));
    }

    static void
    erase_i_j(f_t& a, std::size_t i, std::size_t j)
    {
      base_array_type b = flex_as_base_array(a);
      if (i >= b.size()) scitbx::boost_python::raise_index_error();
      if (j >= b.size()) scitbx::boost_python::raise_index_error();
      b.erase(&b[i], &b[j]);
      a.resize(flex_grid<>(b.size()));
    }

    static void
    resize_1d_1(f_t& a, std::size_t sz)
    {
      base_array_type b = flex_as_base_array(a);
      b.resize(sz);
      a.resize(flex_grid<>(b.size()));
    }

    static void
    resize_1d_2(f_t& a, std::size_t sz, e_t const& x)
    {
      base_array_type b = flex_as_base_array(a);
      b.resize(sz, x);
      a.resize(flex_grid<>(b.size()));
    }

    static void
    resize_flex_grid_1(f_t& a, flex_grid<> const& grid)
    {
      a.resize(grid);
    }

    static void
    resize_flex_grid_2(f_t& a, flex_grid<> const& grid, e_t const& x)
    {
      a.resize(grid, x);
    }

    static void
    clear(f_t& a)
    {
      base_array_type b = flex_as_base_array(a);
      b.clear();
      a.resize(flex_grid<>(b.size()));
    }

    static void
    append(f_t& a, f_t const& other)
    {
      base_array_type b = flex_as_base_array(a);
      assert_0_based_1d(other.accessor());
      b.insert(b.end(), other.begin(), other.end());
      a.resize(flex_grid<>(b.size()));
    }

    static boost::python::object
    indices(f_t const& a)
    {
      if (!a.check_shared_size()) raise_shared_size_mismatch();
      return scitbx::boost_python::range(a.size());
    }

    static flex_items<e_t>
    items(f_t const& a)
    {
      if (!a.check_shared_size()) raise_shared_size_mismatch();
      return flex_items<e_t>(a);
    }

    static f_t
    select(f_t const& a, flex_bool const& flags)
    {
      if (!a.check_shared_size()) raise_shared_size_mismatch();
      assert_0_based_1d(a.accessor());
      assert_0_based_1d(flags.accessor());
      if (a.size() != flags.size()) {
        raise_incompatible_arrays();
      }
      std::size_t n = 0;
      std::size_t i;
      for(i=0;i<flags.size();i++) if (flags[i]) n++;
      base_array_type result;
      result.reserve(n);
      for(i=0;i<flags.size();i++) if (flags[i]) result.push_back(a[i]);
      return f_t(result, flex_grid<>(result.size()));
    }

    static f_t
    shuffle(f_t const& a, flex_size_t const& permutation)
    {
      if (!a.check_shared_size()) raise_shared_size_mismatch();
      assert_0_based_1d(a.accessor());
      assert_0_based_1d(permutation.accessor());
      if (a.size() != permutation.size()) {
        raise_incompatible_arrays();
      }
      base_array_type result;
      if (a.size()) {
        result.resize(a.size(), a[0]); // avoid requirement that e_t is
        for(std::size_t i=1;i<a.size();i++) {  // default constructible
          std::size_t j = permutation[i];
          SCITBX_ASSERT(j < a.size());
          result[j] = a[i];
        }
      }
      return f_t(result, a.accessor());
    }

    static flex_bool
    invert_a(flex_bool const& a) { return !a; }

    static flex_bool
    and_a_a(flex_bool const& a1, flex_bool const& a2) { return a1 && a2; }

    static flex_bool
    or_a_a(flex_bool const& a1, flex_bool const& a2) { return a1 || a2; }

    static flex_bool
    iand_a_a(flex_bool a1, flex_bool const& a2)
    {
      if (a1.accessor() != a2.accessor()) {
        raise_incompatible_arrays();
      }
      for(std::size_t i=0;i<a1.size();i++) if(!a2[i]) a1[i] = false;
      return a1;
    }

    static flex_bool
    ior_a_a(flex_bool a1, flex_bool const& a2)
    {
      if (a1.accessor() != a2.accessor()) {
        raise_incompatible_arrays();
      }
      for(std::size_t i=0;i<a1.size();i++) if(a2[i]) a1[i] = true;
      return a1;
    }

    static flex_bool
    iand_a_s(flex_bool a1, bool a2)
    {
      if (!a2) std::fill(a1.begin(), a1.end(), false);
      return a1;
    }

    static flex_bool
    ior_a_s(flex_bool a1, bool a2)
    {
      if (a2) std::fill(a1.begin(), a1.end(), true);
      return a1;
    }

    static std::size_t
    count(flex_bool const& a1, bool a2)
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

    static f_t neg_a(f_t const& a) { return -a; }
    static f_t add_a_a(f_t const& a1, f_t const& a2) { return a1 + a2; }
    static f_t sub_a_a(f_t const& a1, f_t const& a2) { return a1 - a2; }
    static f_t mul_a_a(f_t const& a1, f_t const& a2) { return a1 * a2; }
    static f_t div_a_a(f_t const& a1, f_t const& a2) { return a1 / a2; }
    static f_t mod_a_a(f_t const& a1, f_t const& a2) { return a1 % a2; }
    static f_t add_a_s(f_t const& a1, e_t const& a2) { return a1 + a2; }
    static f_t sub_a_s(f_t const& a1, e_t const& a2) { return a1 - a2; }
    static f_t mul_a_s(f_t const& a1, e_t const& a2) { return a1 * a2; }
    static f_t div_a_s(f_t const& a1, e_t const& a2) { return a1 / a2; }
    static f_t mod_a_s(f_t const& a1, e_t const& a2) { return a1 % a2; }
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
    static e_t sum_a(f_t const& a) { return sum(a); }
    static e_t sum_sq_a(f_t const& a) { return sum_sq(a); }
    static e_t product_a(f_t const& a) { return product(a); }
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

    static versa<double, flex_grid<> >
    real_complex(versa<std::complex<double>, flex_grid<> > const& a)
    {
      return real(a);
    }

    static versa<double, flex_grid<> >
    imag_complex(versa<std::complex<double>, flex_grid<> > const& a)
    {
      return imag(a);
    }

    static versa<double, flex_grid<> >
    abs_complex(versa<std::complex<double>, flex_grid<> > const& a)
    {
      shared_plain<double> result(a.size(), init_functor_null<double>());
      for(std::size_t i=0;i<a.size();i++) result[i] = std::abs(a[i]);
      return versa<double, flex_grid<> >(result, a.accessor());
    }

    static versa<double, flex_grid<> >
    arg_complex_2(versa<std::complex<double>, flex_grid<> > const& a, bool deg)
    {
      shared_plain<double> result(a.size(), init_functor_null<double>());
      for(std::size_t i=0;i<a.size();i++) {
        result[i] = std::arg(a[i]);
        if (deg) result[i] /= constants::pi_180;
      }
      return versa<double, flex_grid<> >(result, a.accessor());
    }

    static versa<double, flex_grid<> >
    arg_complex_1(versa<std::complex<double>, flex_grid<> > const& a)
    {
      return arg_complex_2(a, false);
    }

    static versa<double, flex_grid<> >
    norm_complex(versa<std::complex<double>, flex_grid<> > const& a)
    {
      return norm(a);
    }

    static versa<std::complex<double>, flex_grid<> >
    polar_complex_3(
      versa<double, flex_grid<> > const& rho,
      versa<double, flex_grid<> > const& theta,
      bool deg)
    {
      if (rho.accessor() != theta.accessor()) {
        raise_incompatible_arrays();
      }
      shared_plain<std::complex<double> > result(
        rho.size(), init_functor_null<std::complex<double> >());
      if (deg) {
        for(std::size_t i=0;i<rho.size();i++) {
          result[i] = std::polar(rho[i], theta[i] * constants::pi_180);
        }
      }
      else {
        for(std::size_t i=0;i<rho.size();i++) {
          result[i] = std::polar(rho[i], theta[i]);
        }
      }
      return versa<std::complex<double>, flex_grid<> >(result, rho.accessor());
    }

    static versa<std::complex<double>, flex_grid<> >
    polar_complex_2(
      versa<double, flex_grid<> > const& rho,
      versa<double, flex_grid<> > const& theta)
    {
      return polar_complex_3(rho, theta, false);
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

      class_<flex_items<ElementType> >((python_name+"_items").c_str())
        .def(init<>())
        .def(init<f_t const&>())
        .def("__len__", &flex_items<ElementType>::size)
        .def("__getitem__", &flex_items<ElementType>::getitem)
      ;

      return class_f_t(python_name.c_str())
        .def(init<>())
        .def(init<flex_grid<> const&, optional<ElementType const&> >())
        .def(init<std::size_t, optional<ElementType const&> >())
        .def(init<shared_plain<ElementType> const&>())
        .def("accessor", accessor)
        .def("nd", nd)
        .def("origin", origin)
        .def("grid", grid)
        .def("last", last_0)
        .def("last", last_1)
        .def("layout", layout)
        .def("is_0_based", is_0_based)
        .def("is_padded", is_padded)
        .def("id", id)
        .def("size", size)
        .def("__len__", size)
        .def("capacity", capacity)
        .def("__getitem__", getitem_1d, GetitemReturnValuePolicy())
        .def("__getitem__", getitem_flex_grid, GetitemReturnValuePolicy())
        .def("__setitem__", setitem_1d)
        .def("__setitem__", setitem_flex_grid)
        .def("front", front, GetitemReturnValuePolicy())
        .def("back", back, GetitemReturnValuePolicy())
        .def("fill", fill)
        .def("deep_copy", deep_copy)
        .def("shallow_copy", shallow_copy)
        .def("as_1d", as_1d)
        .def("assign", assign)
        .def("push_back", push_back)
        .def("pop_back", pop_back)
        .def("append", push_back)
        .def("insert", insert_i_x)
        .def("insert", insert_i_n_x)
        .def("erase", erase_i)
        .def("erase", erase_i_j)
        .def("resize", resize_1d_1)
        .def("resize", resize_1d_2)
        .def("resize", resize_flex_grid_1)
        .def("resize", resize_flex_grid_2)
        .def("clear", clear)
        .def("append", append)
        .def("indices", indices)
        .def("items", items)
        .def("select", select)
        .def("shuffle", shuffle)
      ;
    }

    static class_f_t
    ordered(std::string const& python_name,
            boost::python::object const& flex_root_scope)
    {
      using namespace boost::python;
      using boost::python::def; // works around gcc 2.96 bug
      {
        scope local_scope(flex_root_scope);
        def("order", order_a_a);
      }
      return plain(python_name);
    }

    static class_f_t
    logical(std::string const& python_name,
            boost::python::object const& flex_root_scope)
    {
      return ordered(python_name, flex_root_scope)
        .def("__invert__", invert_a)
        .def("__and__", and_a_a)
        .def("__or__", or_a_a)
        .def("__iand__", iand_a_a)
        .def("__ior__", ior_a_a)
        .def("__iand__", iand_a_s)
        .def("__ior__", ior_a_s)
        .def("count", count)
      ;
    }

    static class_f_t
    numeric_common(std::string const& python_name,
                   boost::python::object const& flex_root_scope)
    {
      using namespace boost::python;
      using boost::python::def; // works around gcc 2.96 bug
      {
        scope local_scope(flex_root_scope);
        def("sum", sum_a);
        def("sum_sq", sum_sq_a);
        def("product", product_a);
      }
      return plain(python_name)
        .def("__neg__", neg_a)
        .def("__add__", add_a_a)
        .def("__sub__", sub_a_a)
        .def("__mul__", mul_a_a)
        .def("__div__", div_a_a)
        .def("__add__", add_a_s)
        .def("__sub__", sub_a_s)
        .def("__mul__", mul_a_s)
        .def("__div__", div_a_s)
        .def("__iadd__", iadd_a_s)
        .def("__isub__", isub_a_s)
        .def("__imul__", imul_a_s)
        .def("__idiv__", idiv_a_s)
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
      using namespace boost::python;
      using boost::python::def; // works around gcc 2.96 bug
      {
        scope local_scope(flex_root_scope);
        def("order", order_a_a);
        def("abs", abs_a);
        def("pow2", pow2_a);
        def("min_index", min_index_a);
        def("max_index", max_index_a);
        def("min", min_a);
        def("max", max_a);
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
      ;
    }

    static class_f_t
    numeric(std::string const& python_name,
            boost::python::object const& flex_root_scope)
    {
      using namespace boost::python;
      using boost::python::def; // works around gcc 2.96 bug
      {
        scope local_scope(flex_root_scope);
        def("pow", pow_a_s);
        def("fmod", fmod_a_s);
        def("atan2", atan2_a_a);
        def("acos", acos_a);
        def("cos", cos_a);
        def("tan", tan_a);
        def("asin", asin_a);
        def("cosh", cosh_a);
        def("tanh", tanh_a);
        def("atan", atan_a);
        def("exp", exp_a);
        def("sin", sin_a);
        def("fabs", fabs_a);
        def("log", log_a);
        def("sinh", sinh_a);
        def("ceil", ceil_a);
        def("floor", floor_a);
        def("log10", log10_a);
        def("sqrt", sqrt_a);
        def("mean", mean_a);
        def("mean_sq", mean_sq_a);
        def("mean_weighted", mean_weighted_a_a);
        def("mean_sq_weighted", mean_sq_weighted_a_a);
      }
      return numeric_no_pow(python_name, flex_root_scope);
    }

    static class_f_t
    integer(std::string const& python_name,
            boost::python::object const& flex_root_scope)
    {
      return numeric_no_pow(python_name, flex_root_scope)
        .def("__mod__", mod_a_a)
        .def("__mod__", mod_a_s)
        .def("__imod__", imod_a_s)
      ;
    }

    static class_f_t
    complex(std::string const& python_name,
            boost::python::object const& flex_root_scope)
    {
      using namespace boost::python;
      using boost::python::def; // works around gcc 2.96 bug
      {
        scope local_scope(flex_root_scope);
        def("real", real_complex);
        def("imag", imag_complex);
        def("abs", abs_complex);
        def("arg", arg_complex_2);
        def("arg", arg_complex_1);
        def("norm", norm_complex);
        def("polar", polar_complex_3);
        def("polar", polar_complex_2);
      }
      return numeric_common(python_name, flex_root_scope);
    }
  };

}}} // namespace scitbx::af::boost_python

#endif // SCITBX_ARRAY_FAMILY_BOOST_PYTHON_FLEX_WRAPPER_H
