#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/list.hpp>
#include <boost/python/enum.hpp>

#include <cctbx/sgtbx/space_group.h>
#include <rstbx/dps_core/dps_core.h>
#include <rstbx/dps_core/direction.h>
#include <rstbx/diffraction/ewald_sphere.h>
#include <rstbx/diffraction/partial_spot_position_partial_H.h>

#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/sort.h>
#include <scitbx/array_family/selections.h>
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
  flex_Direction finegrained_angles;
  scitbx::af::shared<int> n_entries_finegrained;

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

  // Function to construct a more fine-grained hemisphere grid about certain psi/phi angles
  // Supply sampling value and also the old sampling value with which the psi/phi angles may have been generated
  // SST_filter_angles is the flex array of supplied psi/phi angles about which to do finegraining
  flex_Direction construct_hemisphere_grid_finegrained(const double& sampling, const double& old_sampling, flex_Direction SST_filter_angles) {
    SCITBX_ASSERT(sampling <= old_sampling);
    finegrained_angles = flex_Direction(); // Stuff to be returned
    int n_filters=SST_filter_angles.size();
    int psi_index_range = int (0.5 + 2.0*old_sampling/sampling);
    double adjusted_psi_incr = 2*old_sampling/psi_index_range;
    //scitbx::af::shared<int> n_entries_finegrained;
    int n_filter_counter=0;
    n_entries_finegrained.push_back(n_filter_counter);
    for (int ii=0; ii<n_filters; ++ii) {
      Direction dir = SST_filter_angles[ii];
      double max_psi=dir.psi+old_sampling;
      double min_psi=dir.psi-old_sampling;
      double max_phi=dir.phi+old_sampling;
      double min_phi=dir.phi-old_sampling;
      for (int x=0; x <=psi_index_range; ++x) {
        double psi=min_psi + x*adjusted_psi_incr;
        if (psi > scitbx::constants::pi) {
          double eps=1E-4; psi=scitbx::constants::pi-eps;
        }
        if (psi==0) {
          double phi=0;
          finegrained_angles.push_back(Direction(psi,phi));
          n_filter_counter++;
        }
        else {
          int phi_index_range = int (0.5 + 2.0*old_sampling/sampling) ;
          double adjusted_phi_incr = 2*old_sampling/phi_index_range;
          for (int y=0; y <= phi_index_range; ++y) {
            double phi = min_phi + y*adjusted_phi_incr;
            finegrained_angles.push_back(Direction(psi,phi));
            n_filter_counter++;
          }

        }
      }
      n_entries_finegrained.push_back(n_filter_counter);
    }
    return finegrained_angles;
  }

  // coarse grid search ported to c++ for speedup 
  boost::python::tuple coarse_grid_search_cpp(flex_Direction SST_angles, 
                          scitbx::af::shared <double> unique_cell_dimensions, 
                          scitbx::af::const_ref<scitbx::vec3<double> > reciprocal_lattice_vectors) {
    //std::cout << "Wow I am in coarse grid search" << std::endl;
    scitbx::af::shared <scitbx::vec3 < double> > vectors;
    scitbx::af::shared <double> function_values;
    scitbx::af::shared <Direction> SST_all_angles;
    for (int i=0; i < SST_angles.size(); ++i) {
      for (int j=0; j < unique_cell_dimensions.size(); ++j)  {
        scitbx::vec3<double> v = scitbx::vec3<double>(SST_angles[i].dvec[0]*unique_cell_dimensions[j],
                                  SST_angles[i].dvec[1]*unique_cell_dimensions[j],
                                  SST_angles[i].dvec[2]*unique_cell_dimensions[j]);

        double multiplier = 2*scitbx::constants::pi;
        scitbx::af::shared <double> two_pi_S_dot_v = dot_a_s(reciprocal_lattice_vectors,v,multiplier);
        //scitbx::af::shared <double > cosines = all_cos(two_pi_S_dot_v);
        //double f= sum(cosines);
        double f = sum_of_all_cosines(two_pi_S_dot_v);
        vectors.push_back(v);
        function_values.push_back(f);
        SST_all_angles.push_back(SST_angles[i]);
        //std::cout << SST_angles[i].dvec[0]  <<" " << v[0] <<" "<< two_pi_S_dot_v[0]<<" "<<unique_cell_dimensions[j] << std::endl;
        //std::cout << SST_angles[i].dvec[0] << SST_angles[i].dvec[1] << SST_angles[i].dvec[2] << unique_cell_dimensions[j] << v[0]<<v[1]<<v[2]<<std::endl;
      }
    }
    return boost::python::make_tuple(vectors, function_values, SST_all_angles);
  } // coarse_grid_search end

  boost::python::tuple fine_grid_search_cpp(flex_Direction SST_finegrained_angles,
                            scitbx::af::shared <double> unique_cell_dimensions,
                            scitbx::af::const_ref <scitbx::vec3<double> > reciprocal_lattice_vectors) {

    //std::cout << "Wow I am in fine grid search" << std::endl;
    scitbx::af::shared <scitbx::vec3 < double> > vectors;
    scitbx::af::shared <double> function_values;
    int top_n_values = 1; // Number of top scoring vectors to return in each coarse grid per unique dim
    for (int i=0; i<n_entries_finegrained.size()-1; ++i) {
      int start = n_entries_finegrained[i];
      int end = n_entries_finegrained[i+1];
      for (int j=0; j < unique_cell_dimensions.size(); ++j) {
        scitbx::af::shared<scitbx::vec3<double> > tmp_vectors;
        scitbx::af::shared<double> tmp_function_values;
        for (int k=start; k<end; ++k) {
          scitbx::vec3<double> v = scitbx::vec3<double>(SST_finegrained_angles[k].dvec[0]*unique_cell_dimensions[j],
                                  SST_finegrained_angles[k].dvec[1]*unique_cell_dimensions[j],
                                  SST_finegrained_angles[k].dvec[2]*unique_cell_dimensions[j]);
          double multiplier = 2*scitbx::constants::pi;
          scitbx::af::shared <double> two_pi_S_dot_v = dot_a_s(reciprocal_lattice_vectors,v,multiplier);
          //scitbx::af::shared <double > cosines = all_cos(two_pi_S_dot_v);
          //double f= sum(cosines); 
          double f = sum_of_all_cosines(two_pi_S_dot_v);
          tmp_vectors.push_back(v);
          tmp_function_values.push_back(f);
        }
        // Need to sort things here 
        scitbx::af::shared <std::size_t> perm=scitbx::af::sort_permutation(tmp_function_values.const_ref(), true);
        scitbx::af::const_ref<std::size_t> p = perm.const_ref();
        tmp_vectors = scitbx::af::select(tmp_vectors.const_ref(), p);
        tmp_function_values = scitbx::af::select(tmp_function_values.const_ref(), p);
        for (int z=0; z<top_n_values; ++z) {
          vectors.push_back(tmp_vectors[z]);
          function_values.push_back(tmp_function_values[z]);
        } //z
      }// unique cell 
    } // finegrained n_entries
    return boost::python::make_tuple(vectors, function_values);

  } // fine_grid_search end


  // Copying this function from flex_vec3_double.cpp because not sure how to take dot product in c++ for flex_vec3_doubles
  scitbx::af::shared<double>
  dot_a_s(
    scitbx::af::const_ref<scitbx::vec3<double> > const& lhs,
    scitbx::vec3<double> rhs,
    double multiplier)
  {
    scitbx::af::shared<double> result((scitbx::af::reserve(lhs.size())));
    for(std::size_t i=0;i<lhs.size();i++) {
      result.push_back(lhs[i] * rhs*multiplier);
    }
    return result;
  }

  // Copying this function since not sure how to take cosines for flex in c++
  scitbx::af::shared <double> all_cos(scitbx::af::shared<double> mylist)
  {
    scitbx::af::shared<double> result((scitbx::af::reserve(mylist.size())));
    for (std::size_t i=0;  i<mylist.size(); ++i) {
      result.push_back(std::cos(mylist[i]));
      }
      return result;
    }

  // Sum over all elements in a list
  double sum(scitbx::af::shared<double> mylist)
  {
    double result = 0.0;
    for (std::size_t i=0;  i<mylist.size(); ++i) {
      result +=mylist[i];
      }
      return result;
    }

  // Sum over values after taking their cosines
  double sum_of_all_cosines(scitbx::af::shared<double> mylist) {
    double result = 0.0;
    for (std::size_t i=0; i<mylist.size(); ++i) {
      result += std::cos(mylist[i]);
    }
    return result;
  }

};

static boost::python::tuple
observed_indices_and_angles_from_rotation_angles_range(rotation_angles& ra,
  double const& phi_start_rad,double const& phi_end_rad,
  const scitbx::af::shared<cctbx::miller::index<> >& indices){
    // This is going to require some revision to assure it works in an arbitrary
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
    scitbx::af::shared<double> return_angle_full_width_rad;
    scitbx::af::shared<double> return_lorentz_factor;
    scitbx::af::shared<scitbx::vec3<double> > return_s;

    if (rp.use_gh1982a){
      for (int ihkl = 0; ihkl < observed_indices.size(); ++ihkl) {
        if (rp( observed_indices[ihkl], observed_angles[ihkl] )) {
           scitbx::vec2<double> xy = rp.get_prediction();
           scitbx::vec3<double> s = rp.get_s();
           return_indices.push_back(observed_indices[ihkl]);
           return_angle_rad.push_back(observed_angles[ihkl]);
           return_fast_px.push_back(xy[0]);
           return_slow_px.push_back(xy[1]);
           return_lorentz_factor.push_back(rp.lorentz_factor());
           return_angle_full_width_rad.push_back(rp.get_full_width());
           return_s.push_back(s);
         }
       }
       return make_tuple(return_indices,return_fast_px,return_slow_px,return_angle_rad,
                         return_lorentz_factor, return_angle_full_width_rad, return_s);
    }
    for (int ihkl = 0; ihkl < observed_indices.size(); ++ihkl) {
       if (rp( observed_indices[ihkl], observed_angles[ihkl] )) {
         scitbx::vec2<double> xy = rp.get_prediction();
         scitbx::vec3<double> s = rp.get_s();

                 return_indices.push_back(observed_indices[ihkl]);
                 return_angle_rad.push_back(observed_angles[ihkl]);
                 return_fast_px.push_back(xy[0]);
                 return_slow_px.push_back(xy[1]);
           return_s.push_back(s);
       }
     }
      return make_tuple(return_indices,return_fast_px,return_slow_px,return_angle_rad, return_s);
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
      .def("getXyzData",&dps_core::getXyzData)
      .def("fft_result",fft_result)
      .def("setSolutions",&dps_core::setSolutions)
      .def("set_presorted_solutions",&dps_core::set_presorted_solutions)
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
      .def_readonly("granularity",&dps_core::granularity)
      .def_readonly("amax",&dps_core::amax)

    ;

    class_<Direction>("Direction", init<const point &>())
      .def(init<const double &, const double &>((arg("psi"),arg("phi"))))
      .def_readonly("kmax",&Direction::kmax)
      .add_property("kval",make_getter(&Direction::kval, rbv()),
                           make_setter(&Direction::kval, dcp()))
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
    class_<Directional_FFT>("Directional_FFT", no_init)
      .def(init<const Direction&, const af::shared<scitbx::vec3<double> >&,
                const double&, const double&,
                const sztype&>((arg("angle"),arg("xyzdata"),
                                arg("granularity"),arg("amax"),
                                arg("F0_cutoff"))))
      .def_readonly("pmin",&Directional_FFT::pmin)
      .def_readonly("delta_p",&Directional_FFT::delta_p)
      .def_readonly("fft_result",&Directional_FFT::fft_result)
      .def("kval",&Directional_FFT::kval)
      .def("kmax",&Directional_FFT::kmax)
   ;

    class_<SimpleSamplerTool >("SimpleSamplerTool", init<const double &>())
      .add_property("incr",make_getter(&SimpleSamplerTool::incr, rbv()),
                           make_setter(&SimpleSamplerTool::incr, dcp()))
      .add_property("angles",make_getter(&SimpleSamplerTool::angles, rbv()),
                             make_setter(&SimpleSamplerTool::angles, dcp()))
      .add_property("finegrained_angles",make_getter(&SimpleSamplerTool::finegrained_angles, rbv()),
                             make_setter(&SimpleSamplerTool::finegrained_angles, dcp()))
      .add_property("n_entries_finegrained",make_getter(&SimpleSamplerTool::n_entries_finegrained, rbv()),
                             make_setter(&SimpleSamplerTool::n_entries_finegrained, dcp()))
      .def("construct_hemisphere_grid",&SimpleSamplerTool::construct_hemisphere_grid)
      .def("construct_hemisphere_grid_finegrained",&SimpleSamplerTool::construct_hemisphere_grid_finegrained,
                                                   "Function to construct a more fine-grained hemisphere grid about certain psi/phi angles\n"
                                                   "@params sampling : grid value to construct finegrained grid with\n"
                                                   "@params old_sampling : old grid value that was used for coarse-grained grid\n"
                                                   "@params SST_filter_angles: flex array of supplied psi/phi angles (flex.Direction) about which to do finegraining\n"
      )
      .def ("coarse_grid_search_cpp", &SimpleSamplerTool::coarse_grid_search_cpp)
      .def ("fine_grid_search_cpp", &SimpleSamplerTool::fine_grid_search_cpp)
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
                                  const reflection_prediction::sensor_type &>
      ((arg("axis"), arg("s0"), arg("ub"), arg("sensor"))))
      .def("__call__", & reflection_prediction::operator())
      .def("get_prediction", & reflection_prediction::get_prediction)
      .def("get_s", & reflection_prediction::get_s)
      .def("set_rocking_curve", & reflection_prediction::set_rocking_curve)
      .def("set_mosaicity", & reflection_prediction::set_mosaicity,
            (arg("mos"),arg("degrees")))
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

      .def("vec3", &scattering_list::mm_coord) // when used as a spot positions container
      .def("hkl", &scattering_list::reflections) // when used as a spot positions container
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
