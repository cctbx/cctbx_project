#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>

#include <rstbx/dps_core/dps_core.h>
#include <rstbx/dps_core/direction.h>

using namespace boost::python;

namespace rstbx { namespace boost_python { namespace {

  boost::python::tuple
  foo()
  {
    return boost::python::make_tuple(1,2,3,4);
  }

  Direction
  fft_result(dps_core& ai,Direction& angle){
      fftptr dfft( ai.fft_factory(angle) );
      dfft->extract_directional_properties(true);
      return angle;
  }

  void
  init_module() {
    using namespace boost::python;

    typedef return_value_policy<return_by_value> rbv;
    typedef default_call_policies dcp;

    def("foo", &foo);

    class_<dps_core>("dps_core",init<>())
      .def("setMaxcell",&dps_core::setMaxcell)
      .def("setXyzData",&dps_core::setXyzData)
      .def("getXyzSize",&dps_core::getXyzSize)
      .def("fft_result",fft_result)
      .def("setSolutions",&dps_core::setSolutions)
      .def("getSolutions",&dps_core::getSolutions)
      .def("n_candidates",&dps_core::n_candidates)
      .def("__getitem__",&dps_core::candidate)
      .def("setOrientation",&dps_core::setOrientation)
      .def("set_orientation_direct_matrix",
          &dps_core::set_orientation_direct_matrix)
      .def("set_orientation_reciprocal_matrix",
          &dps_core::set_orientation_reciprocal_matrix)
      .def("getOrientation",&dps_core::getOrientation)
      .def("rmsdev",&dps_core::rmsdev)
      .def("hklobserved",(hkllistmm(dps_core::*)()const) &dps_core::hklobserved)
      .def("hklobserved",
(hkllistmm (dps_core::*)(const pointlistmm&)const) &dps_core::hklobserved)
      .def("observed",&dps_core::observed)
    ;

    class_<Direction>("Direction", init<const point &>())
      .def(init<const double &, const double &>((arg_("psi"),arg_("phi"))))
      .def_readonly("kmax",&Direction::kmax)
      .def_readonly("kval",&Direction::kval)
      .def_readonly("kval0",&Direction::kval0)
      .def_readonly("kval2",&Direction::kval2)
      .def_readonly("kval3",&Direction::kval3)
      .def_readonly("psi",&Direction::psi)
      .def_readonly("phi",&Direction::phi)
      .def_readonly("m",&Direction::m)
      .def_readonly("delta_p",&Direction::delta_p)
      .def("getff",&Direction::getff)
      .add_property("dvec",make_getter(&Direction::dvec, rbv()))
      .add_property("real",make_getter(&Direction::uc_length, rbv()))
      .def("is_nearly_collinear",&Direction::is_nearly_collinear)
   ;

  }

}}} // namespace omptbx::boost_python::<anonymous>

BOOST_PYTHON_MODULE(rstbx_ext)
{
  rstbx::boost_python::init_module();
}
