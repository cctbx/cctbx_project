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
#include <vector>
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

struct backplane {
 public:
  int boxnbg;
  double boxmean, boxvar, boxstd;
  double Sum_x,Sum_x2;
  backplane(){
  //    cout<<"base constructor"<<endl;
clear(); }
  virtual void accumulate (const int&, const int&, const int& px){
    //cout<<"base accumulate"<<endl;
    Sum_x += px;
    Sum_x2 += (double)(px)*px;
    boxnbg++;
  }
  virtual void clear(){
    boxnbg = 0;
    boxmean = boxvar = boxstd = 0.;
    Sum_x = Sum_x2 = 0.;
  }
  virtual void reinitialize(const int&, const int&){
    //cout<<"base reinitialize"<<endl;
    clear();
  }
  virtual void finish(){
    //cout<<"base finish"<<endl;
    boxmean = Sum_x/boxnbg;
    boxvar = Sum_x2/boxnbg - boxmean*boxmean;
    boxstd = std::sqrt(boxvar);
  }
  virtual inline double localmean (const int&, const int&) {
      //cout<<"base mean "<<boxmean<<endl;
    return boxmean; }

  // Trivial calculation of detector gain; re-check derivation before commenting in
  //inline double gain () const { //Suggested by Mosflm manual: digitized value = GAIN * Equiv # of photons
  //  double sample_average = Sum_x/boxnbg;
  //  double sample_variance = (1./boxnbg)*(Sum_x2 - sample_average*sample_average);
  //  return std::sqrt(sample_variance) / sample_average; // is mean_squared / mean
  //}

  virtual ~backplane(){}
};

struct corrected_backplane: public backplane {
 private:
  int Sum_p2,Sum_pq,Sum_p,Sum_q2,Sum_q;
  double Sum_xp,Sum_xq;
  int xstart,ystart;
  double a,b,c;
  std::vector<int> rho_cache;
  std::vector<int> p_cache;
  std::vector<int> q_cache;
  double rmsd;
  double p,q; //temporary values
 public:
  corrected_backplane(const int& xst, const int& yst):
    xstart(xst),ystart(yst) {
    //cout<<"corrected constructor"<<endl;
    clear();
  }
  void reinitialize(const int& xst, const int& yst){
    //cout<<"corrected reinitialize"<<endl;
    xstart = xst; ystart = yst; clear();
  }
  inline void clear(){
    backplane::clear();
    Sum_p2=0;Sum_pq=0;Sum_p=0;Sum_q2=0;Sum_q=0;Sum_xp=0;Sum_xq=0;
    rho_cache.clear();
    p_cache.clear();
    q_cache.clear();
    rmsd=0.;
  }
  inline void accumulate(const int& x, const int& y, const int& px){
    //cout<<"corrected accumulate"<<endl;
    backplane::accumulate(x,y,px);
    int p = x-xstart;
    int q = y-ystart;
    Sum_p2+=p*p;
    Sum_pq+=p*q;
    Sum_p+=p;
    Sum_q2+=q*q;
    Sum_q+=q;
    Sum_xp+=px*p;
    Sum_xq+=px*q;
    rho_cache.push_back(px);p_cache.push_back(p);q_cache.push_back(q);
  }
  inline void finish(){
    //cout<<"corrected finish"<<endl;
    scitbx::mat3<double> rossmann(Sum_p2,Sum_pq,Sum_p,
                          Sum_pq,Sum_q2,Sum_q,
                          Sum_p,Sum_q,boxnbg);
    scitbx::vec3<double> obs(Sum_xp,Sum_xq,Sum_x);
    scitbx::mat3<double> rinv = rossmann.inverse();
    //scitbx::vec3<double> abc = rossmann.inverse()*obs;
    //a=abc[0]; b= abc[1]; c=abc[2];
    a = rinv[0]*Sum_xp + rinv[1]*Sum_xq +rinv[2]*Sum_x;
    b = rinv[3]*Sum_xp + rinv[4]*Sum_xq +rinv[5]*Sum_x;
    c = rinv[6]*Sum_xp + rinv[7]*Sum_xq +rinv[8]*Sum_x;
    for (int v=0; v<boxnbg; ++v){
      double bgobs_bgplane = rho_cache[v] - a*p_cache[v] - b*q_cache[v] -c;
      rmsd +=  bgobs_bgplane*bgobs_bgplane;
    }
    rmsd = std::sqrt(rmsd/boxnbg); //box standard deviation
    boxstd = rmsd;
  }
  inline double localmean(const int&x, const int&y){
    return a*(x-xstart)+b*(y-ystart)+c;
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
