#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <scitbx/histogram.h>

#include <boost/python/tuple.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/overloads.hpp>
#include <scitbx/boost_python/is_polymorphic_workaround.h>

#include <scitbx/wigner/boost_python/wigner3j.h>

namespace scitbx { namespace wigner {

namespace{
  //fast wigner3j
  struct wigner3j_fast_wrapper
  {
    typedef wigner3j_fast < double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("wigner3j_fast", no_init)
        .def( init< int const&
                  >
             (( arg_("max")
             ))
            )
        .def("compute", &w_t::compute)
      ;
    }

  };

  // general wigner3j
  struct wigner3j_wrapper
  {
    typedef wigner3j < double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("wigner3j", no_init)
        .def( init< int const&,
                    int const&,
                    int const&,
                    int const&,
                    int const&,
                    int const&
                  >
             (( arg_("j1"),
                arg_("j2"),
                arg_("j3"),
                arg_("m1"),
                arg_("m2"),
                arg_("m3")
             ))
            )
        .def("check", &w_t::check)
        .def("get_value", &w_t::get_value)
      ;
    }

  };

// with all m's = 0
  struct wigner3j_zero_wrapper
  {
    typedef wigner3j_zero < double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("wigner3j_zero", no_init)
        .def( init< int const&,
                    int const&,
                    int const&
                  >
             (( arg_("j1"),
                arg_("j2"),
                arg_("j3")
             ))
            )
        .def("get_value", &w_t::get_value)
      ;
    }

  };


} // namespace <anonymous>

namespace boost_python{

  void wrap_wigner3j_fast()
  {
   scitbx::wigner::wigner3j_fast_wrapper::wrap();
  }

  void wrap_wigner3j()
  {
   scitbx::wigner::wigner3j_wrapper::wrap();
  }

  void wrap_wigner3j_zero()
  {
   scitbx::wigner::wigner3j_zero_wrapper::wrap();
  }

}

}}
