#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>

#include <scitbx/math/zernike.h>
#include <scitbx/math/2d_zernike_moments.h>

namespace scitbx { namespace math {
namespace {

  struct two_d_voxel_wrapper
  {
    typedef voxel_2d <double> w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("two_d_voxel", no_init)
        .def( init<
                   int const&,
                   double const&,
                   double const&,
                   double const&,
                   scitbx::af::const_ref< scitbx::vec3<double> >
                  >
             ((
                arg("splat_range"),
                arg("external_rmax"),
                arg("dx"),
                arg("fraction"),
                arg("xyz")
             ))
            )
        .def("rmax", &w_t::rmax)
        .def("np", &w_t::np)
        .def("get_image", &w_t::get_image)
        .def("get_value", &w_t::get_value)
      ;
    }
  };

//
//
  struct two_d_grid_wrapper
  {
    typedef grid_2d  < double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("two_d_grid", no_init)
        .def( init<
                   int const&,
                   int const&
                  >
             ((
                arg("np"),
                arg("n_max")
             ))
            )
        .def("get_ss", &w_t::get_ss)
        .def("get_ss", &w_t::get_all_ss)
        .def("clean_space", &w_t::clean_space)
        .def("construct_space_sum",&w_t::construct_space_sum)
      ;
    }
  };

//
//
  struct two_d_moments_wrapper
  {
    typedef zernike_2d_moments  < double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("two_d_zernike_moments", no_init)
        .def( init<
                   grid_2d< double >,
                   int const&
                  >
             ((
                arg("grid"),
                arg("nmax")
             ))
            )
        .def("moments", &w_t::all_moments)
        .def("nm", &w_t::nm)
        .def("get_moment",&w_t::get_moment)
        .def("calc_moments",&w_t::calc_moments)
        .def("update_ss",&w_t::update_ss)
        .def("print_bnmk",&w_t::print_Bnmk)
        .def("zernike_poly",&w_t::zernike_poly)
        .def("zernike_map",&w_t::zernike_map)
      ;
    }
  };

//
//

} //namespace <anonymous>

namespace boost_python {

  void wrap_2d_zernike_mom()
  {
    two_d_voxel_wrapper::wrap();
    two_d_moments_wrapper::wrap();
    two_d_grid_wrapper::wrap();
  }

}}}
