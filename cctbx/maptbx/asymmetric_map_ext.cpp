#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/def.hpp>

#include <scitbx/array_family/boost_python/c_grid_flex_conversions.h>

#include "cctbx/maptbx/asymmetric_map.h"

namespace cctbx { namespace maptbx { namespace {

void init_module()
{
  typedef asymmetric_map w_t;
  using namespace boost::python;
  typedef return_value_policy<return_by_value> rbv;
  typedef scitbx::af::const_ref<double, scitbx::af::flex_grid<> > cm_t;
  typedef scitbx::af::versa<double, scitbx::af::flex_grid<> > m_t;
  class_<w_t,boost::noncopyable>("asymmetric_map", no_init)
    .def(init< const cctbx::sgtbx::space_group_type &, cm_t >
        ( (arg("space_group_type"), arg("density_map")) ))
    .def(init< const cctbx::sgtbx::space_group_type &, m_t,
        const scitbx::af::int3 & >
        ( (arg("space_group_type"), arg("density_map"), arg("grid_size")) ))
    .def("structure_factors", &w_t::structure_factors, (arg("indices")))
    .def("map_for_fft", &w_t::map_for_fft)
    .def("data", &w_t::data, rbv())
  ;

  scitbx::af::boost_python::c_grid_flex_conversions<double,
    scitbx::af::c_interval_grid<3ul> >();
}

}}} // namespace cctbx::maptbx::anonymous

BOOST_PYTHON_MODULE(cctbx_asymmetric_map_ext)
{
  cctbx::maptbx::init_module();
}
