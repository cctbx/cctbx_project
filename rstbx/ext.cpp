#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/enum.hpp>

#include <cctbx/sgtbx/space_group.h>
#include <rstbx/dps_core/dps_core.h>
#include <rstbx/dps_core/direction.h>
#include <rstbx/diffraction/ewald_sphere.h>
#include <rstbx/diffraction/partial_spot_position_partial_H.h>

#include <scitbx/array_family/flex_types.h>
#include <vector>
#include <rstbx/backplane.h>

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

static boost::python::tuple
observed_indices_and_angles_from_rotation_angles_range(rotation_angles& ra,
  double const& phi_start_rad,double const& phi_end_rad,
  const scitbx::af::shared<cctbx::miller::index<> >& indices){
    //This is going to require some revision to assure it works in an arbitrary
    // principle value region for phi_start_rad and phi_end_rad

    scitbx::af::shared<scitbx::vec3<double> > return_indices;
    scitbx::af::shared<double> return_angles_rad;

    for (int ihkl = 0; ihkl < indices.size(); ++ihkl) {
       scitbx::vec3<double> test_index( // convert integer Miller index to double type
         indices[ihkl][0],indices[ihkl][1],indices[ihkl][2]);
       if (ra( test_index )) {
         scitbx::vec2<Angle> intersection_angles = ra.get_intersection_angles();
         for (int iangle = 0; iangle < 2; ++iangle){
           if (intersection_angles[iangle] >= phi_start_rad &&
               intersection_angles[iangle] <= phi_end_rad){
                 return_indices.push_back(test_index);
                 return_angles_rad.push_back(intersection_angles[iangle]);
           }
         }
       }
     }
     return make_tuple(return_indices,return_angles_rad);
}

static boost::python::tuple
rp_predict(reflection_prediction& rp,
  scitbx::af::shared<scitbx::vec3<double> > const& observed_indices,
  scitbx::af::shared<double> const& observed_angles){

    scitbx::af::shared<scitbx::vec3<double> > return_indices;
    scitbx::af::shared<double> return_fast_px;
    scitbx::af::shared<double> return_slow_px;
    scitbx::af::shared<double> return_angle_rad;

    for (int ihkl = 0; ihkl < observed_indices.size(); ++ihkl) {
       if (rp( observed_indices[ihkl], observed_angles[ihkl] )) {
         scitbx::vec2<double> xy = rp.get_prediction();

                 return_indices.push_back(observed_indices[ihkl]);
                 return_angle_rad.push_back(observed_angles[ihkl]);
                 return_fast_px.push_back(xy[0]);
                 return_slow_px.push_back(xy[1]);
       }
     }
     return make_tuple(return_indices,return_fast_px,return_slow_px,return_angle_rad);
}

static af::shared<cctbx::miller::index<> >
full_sphere_indices(cctbx::uctbx::unit_cell const& uc,
                    double const& resolution_limit,
                    cctbx::sgtbx::space_group const& sg){

  cctbx::miller::index<> maxhkl = uc.max_miller_indices(resolution_limit);

  af::shared<cctbx::miller::index<> > present;

  for (int h = -maxhkl[0]; h <= maxhkl[0]; ++h){
    for (int k = -maxhkl[1]; k <= maxhkl[1]; ++k){
      for (int l = -maxhkl[2]; l <= maxhkl[2]; ++l){

        if (h == 0 && k == 0 && l == 0) { continue; }

        if (uc.d(cctbx::miller::index<>(h, k, l)) < resolution_limit){
                    continue;}

        if (sg.is_sys_absent(cctbx::miller::index<>(h, k, l))){
                    continue;}

        present.push_back(cctbx::miller::index<>(h, k, l));
      }
    }
  }
  return present;
}

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

    class_<corrected_backplane>("corrected_backplane",init<const int&, const int&>())
      .def("accumulate",&corrected_backplane::accumulate)
      .def("finish",&corrected_backplane::finish)
      .def("localmean",&corrected_backplane::localmean)
      .def_readonly("rmsd",&corrected_backplane::boxstd)
    ;

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
      .def("observed_indices_and_angles_from_angle_range",
            &observed_indices_and_angles_from_rotation_angles_range,
           (arg("phi_start_rad"),arg("phi_end_rad"),arg("indices")))
    ;

    class_<reflection_prediction>("reflection_prediction",
                                  init<const scitbx::vec3<double> &,
                                  const scitbx::vec3<double> &,
                                  const scitbx::mat3<double> &,
                                  const scitbx::vec3<double> &,
                                  const scitbx::vec3<double> &,
                                  const scitbx::vec3<double> &,
                                  const double &, const double &,
                                  const double &, const double &>
      ((arg("axis"), arg("s0"), arg("ub"),
        arg("origin"), arg("fast"), arg("slow"),
        arg("f_min"), arg("f_max"),
        arg("s_min"), arg("s_max"))))
      .def("__call__", & reflection_prediction::operator())
      .def("get_prediction", & reflection_prediction::get_prediction)
      .def("predict",
            &rp_predict,
           (arg("observed_indices"),arg("observed_angles")));

    def("full_sphere_indices",&full_sphere_indices,
      (arg("unit_cell"), arg("resolution_limit"),
       arg("space_group")));

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
