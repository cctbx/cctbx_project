#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <mmtbx/rotamer/fit.h>

namespace mmtbx { namespace rotamer {
namespace {

  void init_module()
  {
    using namespace boost::python;

    class_<fit>("fit",
                 init<double,
                      boost::python::list const&,
                      boost::python::list const&,
                      boost::python::list const&,
                      af::const_ref<double, af::c_grid_padded<3> > const&,
                      af::shared<scitbx::vec3<double> >,
                      cctbx::uctbx::unit_cell const&,
                      af::const_ref<std::size_t> const& ,
                      af::const_ref<double> const&,
                      af::const_ref<double> const&,
                      double,
                      int >
                       ((arg("target_value"),
                         arg("axes"),
                         arg("rotatable_points_indices"),
                         arg("angles_array"),
                         arg("density_map"),
                         arg("all_points"),
                         arg("unit_cell"),
                         arg("selection"),
                         arg("sin_table"),
                         arg("cos_table"),
                         arg("step"),
                         arg("n"))))
      .def("result", &fit::result)
    ;

  }

} // namespace <anonymous>
}} // namespace mmtbx::rotamer

BOOST_PYTHON_MODULE(mmtbx_rotamer_fit_ext)
{
  mmtbx::rotamer::init_module();
}
