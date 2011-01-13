#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>

#include <scitbx/math/zernike.h>
#include <scitbx/math/zernike_moments.h>

namespace scitbx { namespace math {
namespace {

//
//
  struct grid_wrapper
  {
    typedef grid  < double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("sphere_grid", no_init)
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
        .def("construct_space_sum_via_list",&w_t::construct_space_sum_via_list)
        .def("construct_space_sum_via_list",&w_t::construct_space_sum_via_list_only)
        .def("unit_sphere_index", &w_t::unit_sphere_index)
        .def("occupied_sites", &w_t::occupied_sites)
      ;
    }
  };

//
//


  struct moments_wrapper
  {
    typedef moments  < double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("zernike_moments", no_init)
        .def( init<
                   grid<double>,
                   int const&
                  >
             ((
                arg("grid"),
                arg("nmax")
             ))
            )
        .def("moments", &w_t::all_moments)
        .def("get_moment",&w_t::get_moment)
        .def("calc_moments",&w_t::calc_moments)
        .def("update_ss",&w_t::update_ss)
        .def("fnn",&w_t::fnn)
        .def("fnl",&w_t::fnl)
        .def("fnnl",&w_t::fnnl)
      ;
    }
  };

//
//
  struct voxel_wrapper
  {
    typedef voxel < double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("sphere_voxel", no_init)
        .def( init<
                   int const&,
                   int const&,
                   bool const&,
                   bool const&,
                   double const&,
                   double const&,
                   double const&,
                   scitbx::af::const_ref< scitbx::vec3<double> >,
                   scitbx::af::const_ref< double >
                  >
             ((
                arg("np"),
                arg("splat_range"),
                arg("uniform"),
                arg("fixed_dx"),
                arg("external_rmax"),
                arg("fraction"),
                arg("dx"),
                arg("xyz"),
                arg("density")
             ))
            )
        .def("value", &w_t::get_value)
        .def("rmax", &w_t::rmax)
        .def("rg", &w_t::rg)
        .def("map", &w_t::map)
        .def("xyz", &w_t::xyz)
        .def("rotate", &w_t::rotate)
        .def("np",&w_t::np)
        .def("occupied_sites", &w_t::occupied_sites)
        .def("status", &w_t::print_status)
        .def("border", &w_t::border)
      ;
    }
  };

} //namespace <anonymous>

namespace boost_python {

  void wrap_zernike_mom()
  {
    voxel_wrapper::wrap();
    moments_wrapper::wrap();
    grid_wrapper::wrap();
  }

}}}
