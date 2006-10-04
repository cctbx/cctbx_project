#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/copy_const_reference.hpp>

#include <scitbx/math/resample.h>

namespace scitbx { namespace math {

namespace {


  struct non_parametric_bootstrap_wrappers
  {
    typedef resample::non_parametric_bootstrap<double> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("non_parametric_bootstrap", no_init)
        .def(init< scitbx::af::const_ref<double> const&,
                  long const&>
             ((arg_("observations"),
               arg_("seed") )))
        .def("draw", &w_t::draw)
        .def("draw_from_random_jack_knifed_sample",
             &w_t::draw_from_random_jack_knifed_sample )
        ;
    }


  };

  //----------------------------------------------------------

  struct non_parametric_bootstrap_as_int_wrappers
  {
    typedef resample::non_parametric_bootstrap_as_int<std::size_t> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("non_parametric_bootstrap_as_int", no_init)
        .def(init< scitbx::af::const_ref<std::size_t> const&,
                  long const&>
             ((arg_("observations"),
               arg_("seed") )))
        .def("draw", &w_t::draw)
        .def("draw_from_random_jack_knifed_sample",
             &w_t::draw_from_random_jack_knifed_sample )
        ;
    }


  };

  //----------------------------------------------------------

  struct smooth_bootstrap_wrappers
  {
    typedef resample::smooth_bootstrap<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("smooth_bootstrap", no_init)
        .def(init< scitbx::af::const_ref<double> const&,
                  long const&>
             ((arg_("observations"),
               arg_("seed") )))
        .def("draw", &w_t::draw)
        .def("draw_from_random_jack_knifed_sample",
             &w_t::draw_from_random_jack_knifed_sample )
        ;
    }


  };

} // namespace <anonymous>

namespace boost_python {

  void wrap_resample()
  {
    non_parametric_bootstrap_wrappers::wrap();
    non_parametric_bootstrap_as_int_wrappers::wrap();
    smooth_bootstrap_wrappers::wrap();
  }

}}} // namespace scitbx::math::boost_python
