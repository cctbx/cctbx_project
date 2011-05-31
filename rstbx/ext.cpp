#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/enum.hpp>

#include <rstbx/dps_core/dps_core.h>
#include <rstbx/dps_core/direction.h>
#include <rstbx/diffraction/ewald_sphere.h>
#include <rstbx/diffraction/partial_spot_position_partial_H.h>

#include <scitbx/array_family/flex_types.h>

using namespace boost::python;

namespace rstbx {

/* SimpleSamplerTool samples the unit sphere; outputting a list of surface
   points that are spaced apart a specified number of radians.
*/

class SimpleSamplerTool {
 typedef scitbx::af::shared<Direction > flex_Direction;

 public:
  double incr;  //initial directional spacing in radians
  flex_Direction angles;

  SimpleSamplerTool(const double& characteristic_grid):
    // The maximum allowable characteristic grid should be about 0.029 radians,
    // corresponding to the grid sampling used in the Rossman DPS paper;
    // approx 7900 directions; 0.03 seconds in C++ code.
    // But the characteristic grid sampling should be reflective of the problem at hand =
    // approximately: the observed resolution limit /
    //                most conservative (largest) cell estimate

    incr(characteristic_grid) {
    //construct_hemisphere_grid(incr);
  }

  flex_Direction construct_hemisphere_grid(const double& sampling){
    // psi is the equivalent of latitude, measured as an angle from the North pole
    // rounding:
    int psi_index_range = int (0.5 + scitbx::constants::pi_2/sampling);
    // adjust for integral number
    double adjusted_psi_incr = scitbx::constants::pi_2/psi_index_range;
    angles = flex_Direction();
    angles.reserve(4*psi_index_range*psi_index_range);

    for (int x = 0; x <= psi_index_range; ++x){
      double psi = x * adjusted_psi_incr;
      if (psi > scitbx::constants::pi){
        double eps = 1E-4; psi=scitbx::constants::pi-eps;
      }

      // phi is the equivalent of longitude
      if (psi==0){
        double phi=0.;
        angles.push_back(Direction(psi,phi));
      } else {
        int phi_index_range = int (0.5 + 2.*scitbx::constants::pi*std::sin(psi)/sampling);
        double adjusted_phi_incr = 2.*scitbx::constants::pi/phi_index_range;
        for (int y =0; y < phi_index_range; ++y) {
          double phi = y * adjusted_phi_incr;
          angles.push_back(Direction(psi,phi));
        }
      }
    }
    return angles;
  }
};

namespace boost_python { namespace {

  boost::python::tuple
  foo()
  {
    return boost::python::make_tuple(1,2,3,4);
  }

  Direction
  fft_result(dps_core& ai,Direction& angle){
      fftptr dfft( ai.fft_factory(angle) );
      angle.extract_directional_properties(dfft,true);
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
      .def(init<const double &, const double &>((arg("psi"),arg("phi"))))
      .def_readonly("kmax",&Direction::kmax)
      .def_readonly("kval",&Direction::kval)
      .def_readonly("kval0",&Direction::kval0)
      .def_readonly("kval2",&Direction::kval2)
      .def_readonly("kval3",&Direction::kval3)
      .def_readonly("psi",&Direction::psi)
      .def_readonly("phi",&Direction::phi)
      .def_readonly("m",&Direction::m)
      .def_readonly("delta_p",&Direction::delta_p)
      .def("bvec",&Direction::bvec)
      .def("getff",&Direction::getff)
      .add_property("dvec",make_getter(&Direction::dvec, rbv()))
      .add_property("real",make_getter(&Direction::uc_length, rbv()))
      .def("is_nearly_collinear",&Direction::is_nearly_collinear)
   ;

    class_<SimpleSamplerTool >("SimpleSamplerTool", init<const double &>())
      .add_property("incr",make_getter(&SimpleSamplerTool::incr, rbv()),
                           make_setter(&SimpleSamplerTool::incr, dcp()))
      .add_property("angles",make_getter(&SimpleSamplerTool::angles, rbv()),
                             make_setter(&SimpleSamplerTool::angles, dcp()))
      .def("construct_hemisphere_grid",&SimpleSamplerTool::construct_hemisphere_grid)
   ;

    class_<ewald_sphere_base_model>("ewald_sphere_base_model",
      init<const double&, const ewald_sphere_base_model::matrix&, const double&,
           const ewald_sphere_base_model::point&>(
           (arg("limiting_resolution"),arg("orientation"),
            arg("wavelength"),arg("axial_direction"))))
      .def("setH",(void(ewald_sphere_base_model::*)
           (const ewald_sphere_base_model::point&)) &ewald_sphere_base_model::setH)
      .def("setH",(void(ewald_sphere_base_model::*)
           (const cctbx::miller::index<>&)) &ewald_sphere_base_model::setH)
      .add_property("H",make_getter(&ewald_sphere_base_model::H, rbv()))
    ;

    class_<rotation_angles, bases<ewald_sphere_base_model> >("rotation_angles",
      init<const double&, const ewald_sphere_base_model::matrix&, const double&,
           const ewald_sphere_base_model::point&>(
           (arg("limiting_resolution"),arg("orientation"),
            arg("wavelength"),arg("axial_direction"))))
      .def(init<const ewald_sphere_base_model&>())
      .def("__call__", &rotation_angles::operator())
      .def("axis", &rotation_angles::axis)
      .def("offsetdot", &rotation_angles::offsetdot)
      .def("get_intersection_angles", &rotation_angles::get_intersection_angles)
    ;

    class_<partial_spot_position_partial_H, bases<rotation_angles> >(
      "partial_spot_position_partial_H",
      init<const double&, const ewald_sphere_base_model::matrix&, const double&,
           const ewald_sphere_base_model::point&>((arg("limiting_resolution"),
             arg("orientation"),
             arg("wavelength"),arg("axial_direction")
      )))
      .def("__call__", &partial_spot_position_partial_H::operator())
      .def("dangle_", &partial_spot_position_partial_H::dangle_)
    ;

    class_<scattering_list >("scattering_list",
      init<scitbx::af::shared<cctbx::miller::index<> >,
                           const cctbx::crystal_orientation&,
                           scitbx::vec3<double>,
                           scitbx::vec2<double>,
                           const double&,
                           const double&>())
      .def("mm_coord", &scattering_list::mm_coord)
      .def("reflections", &scattering_list::reflections)
    ;


  }

}}} // namespace omptbx::boost_python::<anonymous>

BOOST_PYTHON_MODULE(rstbx_ext)
{
  rstbx::boost_python::init_module();

  // Expose SpotClass to Python
  enum_<rstbx::SpotClass>("SpotClass")
    .value("GOOD",rstbx::GOOD)
    .value("OVERLAP",rstbx::OVERLAP)
    .value("SPINDLE",rstbx::SPINDLE)
    .value("ICE",rstbx::ICE)
    .value("OTHERIMAGE",rstbx::OTHERIMAGE)
    .value("FULL_ENTER",rstbx::FULL_ENTER)
    .value("FULL_EXIT",rstbx::FULL_EXIT)
    .value("ENTER1",rstbx::ENTER1)
    .value("ENTER2",rstbx::ENTER2)
    .value("EXIT3",rstbx::EXIT3)
    .value("EXIT4",rstbx::EXIT4)
    .value("NONE",rstbx::NONE)
    .value("OUTLIER",rstbx::OUTLIER)
    .export_values()
    ;
}
