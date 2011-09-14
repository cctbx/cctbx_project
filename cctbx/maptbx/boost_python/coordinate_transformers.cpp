// copyright (c) Jacob N. Smith; leave this here; use at your whim
#include <cctbx/boost_python/flex_fwd.h>
#include <cctbx/maptbx/coordinate_transformers.h>
#include <scitbx/boost_python/is_polymorphic_workaround.h>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

namespace {
  typedef
    cctbx::maptbx::transform<
      cctbx::cartesian<double>,
      cctbx::grid_point<signed long> >
        tf_c_g;
  typedef
    cctbx::maptbx::transform<
      cctbx::grid_point<signed long>,
      cctbx::cartesian<double> >
        tf_g_c;
}
SCITBX_BOOST_IS_POLYMORPHIC_WORKAROUND(tf_c_g)
SCITBX_BOOST_IS_POLYMORPHIC_WORKAROUND(tf_g_c)

namespace cctbx { namespace maptbx { namespace boost_python {

namespace {

struct transform_coordinate_wrappers {

  typedef scitbx::mat3<double>              matrix3;

  typedef cctbx::fractional<double>                fractional;
  typedef cctbx::grid_point<signed long>           grid_point;
  typedef cctbx::cartesian<double>                 cartesian;
  typedef af::tiny<signed long,dimension_3>        extents;

  typedef transform<fractional,fractional>        frac2frac;
  typedef transform<fractional,cartesian>          frac2cart;
  typedef transform<fractional,grid_point>        frac2grid;

  typedef transform<cartesian,fractional>          cart2frac;
  typedef transform<cartesian,cartesian>          cart2cart;
  typedef transform<cartesian,grid_point>          cart2grid;

  typedef transform<grid_point,fractional>        grid2frac;
  typedef transform<grid_point,cartesian>          grid2cart;
  typedef transform<grid_point,grid_point>        grid2grid;

  static void wrap() {

    using namespace boost::python;

    class_<frac2frac>("frac2frac",init<>())
      .def("__call__",&frac2frac::operator(),arg("coordinate"))
      .def("inverse",&frac2frac::inverse);

    class_<frac2cart>("frac2cart",no_init)
      .def(init<matrix3>())
      .def("__call__",&frac2cart::operator(),arg("coordinate"))
      .def("inverse",&frac2cart::inverse);

    class_<frac2grid>("frac2grid",no_init)
      .def(init<extents>())
      .def("__call__",&frac2grid::operator(),arg("coordinate"))
      .def("fractional_transform",&frac2grid::strange_transform,arg("coordinate"))
      .def("floor_transform",&frac2grid::floor_transform,arg("coordinate"))
      .def("inverse",&frac2grid::inverse);

    class_<cart2frac>("cart2frac",no_init)
      .def(init<matrix3>())
      .def("__call__",&cart2frac::operator(),arg("coordinate"))
      .def("inverse",&cart2frac::inverse);

    class_<cart2cart>("cart2cart",init<>())
      .def("__call__",&cart2cart::operator(),arg("coordinate"))
      .def("inverse",&cart2cart::inverse);

    class_<cart2grid>("cart2grid",no_init)
      .def(init<matrix3,extents>())
      .def("__call__",&cart2grid::operator(),arg("coordinate"))
      .def("inverse",&cart2grid::inverse);

    class_<grid2frac>("grid2frac",no_init)
      .def(init<extents>())
      .def("__call__",&grid2frac::operator(),arg("coordinate"))
      .def("inverse",&grid2frac::inverse);

    class_<grid2cart>("grid2cart",no_init)
      .def(init<extents,matrix3>())
      .def("__call__",&grid2cart::operator(),arg("coordinate"))
      .def("inverse",&grid2cart::inverse);

    class_<grid2grid>("grid2grid",init<>())
      .def("__call__",&grid2grid::operator(),arg("coordinate"))
      .def("inverse",&grid2grid::inverse);

  }
};

} // namespace <anoymous>

void wrap_coordinate_transformers () {
  transform_coordinate_wrappers::wrap();
}

}}} // namespace cctbx::maptbx::boost_python
