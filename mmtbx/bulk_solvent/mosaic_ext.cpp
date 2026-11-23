#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <mmtbx/bulk_solvent/mosaic.h>
#include <mmtbx/f_model/f_model.h>
#include <mmtbx/bulk_solvent/bulk_solvent.h>

#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/list.hpp>

namespace mmtbx { namespace mosaic {
namespace {

  void init_module()
  {
    using namespace boost::python;
    using boost::python::arg;
    typedef return_value_policy<return_by_value> rbv;

  class_<alg2_tg<> >(
      "alg2_tg")
      .def(init<
         boost::python::list   const&,
         af::const_ref<double> const&>(
             (arg("F"),
              arg("i_obs"))))
      .def("target",   &alg2_tg<>::target)
      .def("gradient", &alg2_tg<>::gradient)
      .def("update",   &alg2_tg<>::update)
   ;

   def("normal_matrix", normal_matrix, (arg("G")));

   def("alg4_solve_bounded", alg4_solve_bounded, (
     arg("a"),
     arg("b"),
     arg("lower"),
     arg("upper"),
     arg("x_start"),
     arg("max_sweeps")=10000,
     arg("eps")=1.e-10));

   def("alg4", alg4, (
     arg("G"),
     arg("f_obs"),
     arg("base"),
     arg("lower"),
     arg("upper"),
     arg("x_start"),
     arg("n_cycles")=5));

  }

} // namespace <anonymous>
}} // namespace mmtbx::mosaic

BOOST_PYTHON_MODULE(mmtbx_mosaic_ext)
{
  mmtbx::mosaic::init_module();
}
