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
        .def("clean_space", &w_t::clean_space)
        .def("construct_space_sum",&w_t::construct_space_sum)
        .def("unit_sphere", &w_t::unit_sphere)
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
        .def("fnn",&w_t::Fnn)
        .def("fnl",&w_t::Fnl)
        .def("fnnl",&w_t::Fnnl)
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
                   scitbx::af::const_ref< scitbx::vec3<double> >
                  >
             ((
                arg("np"),
                arg("splat_range"),
                arg("uniform"),
                arg("fixed_dx"),
                arg("fraction"),
                arg("dx"),
                arg("xyz")
             ))
            )
        .def("xyz2voxel", &w_t::xyz2voxel)
        .def("value", &w_t::get_value)
        .def("rmax", &w_t::rmax)
        .def("map", &w_t::map)
        .def("xyz", &w_t::xyz)
        .def("np",&w_t::np)
        .def("occupied_sites", &w_t::occupied_sites)
        .def("status", &w_t::print_status)
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
