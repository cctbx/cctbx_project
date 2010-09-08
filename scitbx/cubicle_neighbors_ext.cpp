#include <boost/python.hpp>

#include <scitbx/cubicle_neighbors.hpp>

namespace scitbx { namespace cubicle_neighbors_ext {

  struct cubicle_neighbors_wrappers
  {
    typedef cubicle_neighbors<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("cubicle_neighbors", no_init)
        .def(init<
          af::const_ref<vec3<double> > const&,
          double const&,
          double const&>((
            arg("main_sites_cart"),
            arg("cubicle_edge"),
            arg("epsilon")=1.e-6)))
        .def("neighbors_of", &w_t::neighbors_of, (
          arg("other_sites_cart"),
          arg("distance_cutoff_sq")))
      ;
    }
  };

  void
  init_module()
  {
    cubicle_neighbors_wrappers::wrap();
  }

}} // namespace scitbx::cubicle_neighbors_ext

BOOST_PYTHON_MODULE(scitbx_cubicle_neighbors_ext)
{
  scitbx::cubicle_neighbors_ext::init_module();
}
