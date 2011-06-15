#ifndef SCITBX_STL_VECTOR_WRAPPER_H
#define SCITBX_STL_VECTOR_WRAPPER_H

#include <boost/python/class.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_non_const_reference.hpp>
#include <boost/python/slice.hpp>
#include <scitbx/boost_python/slice.h>
#include <scitbx/boost_python/utils.h>
#include <scitbx/boost_python/container_conversions.h>
#include <scitbx/array_family/boost_python/ref_from_array.h>
#include <scitbx/misc/positive_getitem_index.h>
#include <scitbx/error.h>
#include <vector>

#if defined(__sgi) && !defined(__GNUC__)
// for details see <scitbx/array_family/boost_python/flex_fwd.h>
namespace scitbx { namespace af { namespace boost_python {
  struct stl_vector_fwd
  {
    friend void f(const_ref<std::size_t> const&);
  };
}}} // namespace scitbx::af::boost_python
#endif // defined(__sgi) && !defined(__GNUC__)

namespace scitbx { namespace stl { namespace boost_python {

  using scitbx::positive_getitem_index;

  template <typename ElementType,
            typename GetitemReturnValuePolicy
              = boost::python::return_value_policy<
                  boost::python::copy_non_const_reference> >
  struct vector_wrapper
  {
    typedef std::vector<ElementType> w_t;
    typedef ElementType e_t;

    static e_t&
    getitem_1d(w_t& self, long i)
    {
      return self[positive_getitem_index(i, self.size())];
    }

    static void
    setitem_1d(w_t& self, long i, e_t const& x)
    {
      self[positive_getitem_index(i, self.size())] = x;
    }

    static void
    delitem_1d(w_t& self, long i)
    {
      self.erase(self.begin()+positive_getitem_index(i, self.size()));
    }

    static w_t
    getitem_1d_slice(w_t const& self, boost::python::slice const& slice)
    {
      scitbx::boost_python::adapted_slice a_sl(slice, self.size());
      w_t result;
      result.reserve(a_sl.size);
      for(long i=a_sl.start;i!=a_sl.stop;i+=a_sl.step) {
        result.push_back(self[i]);
      }
      return result;
    }

    static void
    delitem_1d_slice(w_t& self, boost::python::slice const& slice)
    {
      scitbx::boost_python::adapted_slice a_sl(slice, self.size());
      SCITBX_ASSERT(a_sl.step == 1);
      self.erase(self.begin()+a_sl.start, self.begin()+a_sl.stop);
    }

    static void
    insert(w_t& self, long i, e_t const& x)
    {
      self.insert(self.begin()+positive_getitem_index(i, self.size()), x);
    }

    static void
    append(w_t& self, e_t const& x) { self.push_back(x); }

    static void
    extend(w_t& self, w_t const& other)
    {
      self.insert(self.end(), other.begin(), other.end());
    }

    static boost::python::tuple
    getinitargs(w_t const& self)
    {
      return boost::python::make_tuple(boost::python::tuple(self));
    }

    static void
    wrap(std::string const& python_name)
    {
      using namespace boost::python;
      class_<w_t, boost::shared_ptr<w_t> >(python_name.c_str())
        .def(init<w_t const&>())
        .def(init<std::size_t const&, optional<e_t const&> >())
        .def("size", &w_t::size)
        .def("__len__", &w_t::size)
        .def("__getitem__", getitem_1d, GetitemReturnValuePolicy())
        .def("__setitem__", setitem_1d)
        .def("__delitem__", delitem_1d)
        .def("__getitem__", getitem_1d_slice)
        .def("__delitem__", delitem_1d_slice)
        .def("clear", &w_t::clear)
        .def("insert", insert)
        .def("append", append)
        .def("extend", extend)
        .enable_pickling()
        .def("__getinitargs__", getinitargs)
      ;

      scitbx::boost_python::container_conversions::from_python_sequence<
        w_t,
        scitbx::boost_python::container_conversions
          ::variable_capacity_policy>();

      using scitbx::array_family::boost_python::ref_from_array;
      ref_from_array<w_t, af::const_ref<ElementType> >();
      ref_from_array<w_t, af::ref<ElementType> >();
    }
  };

}}} // namespace scitbx::stl::boost_python

#endif // SCITBX_STL_VECTOR_WRAPPER_H
