#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>

#include <scitbx/math/zernike.h>
#include <scitbx/math/rotation.h>

namespace scitbx { namespace math {
namespace boost_python{

  struct dmatrix_wrapper
  {
    typedef dmatrix < double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("dmatrix", no_init)
        .def( init<
                   int const& ,
                   double const&
                  >
             ((
                arg("l_max"),
                arg("beta")
             ))
            )
        .def("djmn", &w_t::djmn)
      ;
    }
  };

  void
  wrap_dmatrix()
  {
    dmatrix_wrapper::wrap();
  }

// correlation
  struct correlation_wrapper
  {
    typedef correlation < double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("correlation", no_init)
        .def( init<
                   scitbx::math::zernike::nlm_array<double> const&,
                   scitbx::math::zernike::nlm_array<double> const&,
                   int const& ,
                   double const&
                  >
             ((
                arg("f_nlm"),
                arg("m_nlm"),
                arg("l_max"),
                arg("beta")
             ))
            )
        .def("calc_correlation", &w_t::calc_correlation)
        .def("mm_coef", &w_t::mm_coef)
        .def("mhm_coef", &w_t::mhm_coef)
        .def("rotate_moving_obj", &w_t::rotate_moving_obj)
        .def("compare_FM", &w_t::compare_fm)
        .def("set_beta", &w_t::set_beta)
      ;
    }
  };

  void
  wrap_correlation()
  {
    correlation_wrapper::wrap();
  }


}
}}
