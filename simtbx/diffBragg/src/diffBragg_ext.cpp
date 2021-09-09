#include <cctbx/boost_python/flex_fwd.h>
#include <boost/python.hpp>
#include <simtbx/diffBragg/src/diffBragg.h>
#include <simtbx/nanoBragg/nanoBragg.h>
#include <iostream>

using namespace boost::python;
namespace simtbx{
namespace diffBragg{
namespace boost_python { namespace {

  void (simtbx::nanoBragg::diffBragg::*add_diffBragg_spots_A)() = &simtbx::nanoBragg::diffBragg::add_diffBragg_spots;
  void (simtbx::nanoBragg::diffBragg::*add_diffBragg_spots_B)(const nanoBragg::af::shared<size_t>&) = &simtbx::nanoBragg::diffBragg::add_diffBragg_spots;
  void (simtbx::nanoBragg::diffBragg::*add_diffBragg_spots_C)(const nanoBragg::af::shared<size_t>&, boost::python::list per_pix_nominal_hkl)
        = &simtbx::nanoBragg::diffBragg::add_diffBragg_spots;

  static void  set_Ndef(simtbx::nanoBragg::diffBragg& diffBragg, boost::python::tuple const& values) {
      diffBragg.set_ncells_def_values(values);
  }

  static boost::python::tuple get_Ndef(simtbx::nanoBragg::diffBragg const& diffBragg) {
      return boost::python::make_tuple(diffBragg.Nd, diffBragg.Ne, diffBragg.Nf);
  }

  static void  set_Nabc_aniso(simtbx::nanoBragg::diffBragg& diffBragg, boost::python::tuple const& values) {
      diffBragg.isotropic_ncells=false;
      diffBragg.set_ncells_values(values);
  }

  static boost::python::tuple get_Nabc_aniso(simtbx::nanoBragg::diffBragg const& diffBragg) {
      return boost::python::make_tuple(diffBragg.Na, diffBragg.Nb, diffBragg.Nc);
  }

  static void  set_Nabc(simtbx::nanoBragg::diffBragg& diffBragg, double const& value) {
      diffBragg.set_value(9, value); //  9 means Ncells parameter
  }

  static double get_Nabc(simtbx::nanoBragg::diffBragg const& diffBragg) {
      return diffBragg.Na;
  }

  boost::python::tuple get_atom_data(simtbx::nanoBragg::diffBragg &diffBragg){
    boost::python::tuple values;
    values = boost::python::make_tuple(0,0);
    return values;
  }
  //TODO override the set_sources function (or xray_beams) property in nanoBragg in order
  // to set the fpfdp accordingly (if Fhkl2 is set)
  static void set_atom_data(simtbx::nanoBragg::diffBragg & diffBragg,
            boost::python::tuple const& atom_XYZBO){
            /// atom XYZBO is tuple of x coord, ycoords, z coords, B-factors, occupancy
            // e.g. 5 parameters per atom
   boost::python::list atom_x = boost::python::extract<boost::python::list>(atom_XYZBO[0]);
   boost::python::list atom_y = boost::python::extract<boost::python::list>(atom_XYZBO[1]);
   boost::python::list atom_z = boost::python::extract<boost::python::list>(atom_XYZBO[2]);
   boost::python::list atom_B = boost::python::extract<boost::python::list>(atom_XYZBO[3]);
   boost::python::list atom_O = boost::python::extract<boost::python::list>(atom_XYZBO[4]);
    int natoms = boost::python::len(atom_x);
    SCITBX_ASSERT( boost::python::len(atom_y)== natoms);
    SCITBX_ASSERT( boost::python::len(atom_z)== natoms);
    SCITBX_ASSERT( boost::python::len(atom_B)== natoms);
    SCITBX_ASSERT( boost::python::len(atom_O)== natoms);
    diffBragg.db_cryst.atom_data.clear();
    for (int i=0; i < natoms; i++){
        diffBragg.db_cryst.atom_data.push_back( boost::python::extract<double>(atom_x[i]));
        diffBragg.db_cryst.atom_data.push_back(  boost::python::extract<double>(atom_y[i]));
        diffBragg.db_cryst.atom_data.push_back(  boost::python::extract<double>(atom_z[i]));
        diffBragg.db_cryst.atom_data.push_back( boost::python::extract<double>(atom_B[i]));
        diffBragg.db_cryst.atom_data.push_back(  boost::python::extract<double>(atom_O[i]));
    }
  }

  boost::python::tuple get_fpfdp_derivs(simtbx::nanoBragg::diffBragg & diffBragg){
        boost::python::tuple vals;
        vals = boost::python::make_tuple(0,0);
        //TODO implement
        return vals;
  }

  static void set_fpfdp_derivs(simtbx::nanoBragg::diffBragg & diffBragg,
                boost::python::tuple const& derivs) {
    boost::python::list fprime_derivs = boost::python::extract<boost::python::list>(derivs[0]);
    boost::python::list fdblprime_derivs = boost::python::extract<boost::python::list>(derivs[1]);
    int num_fprime = boost::python::len(fprime_derivs);
    int num_fdblprime =   boost::python::len(fdblprime_derivs);
    SCITBX_ASSERT(num_fprime % diffBragg.sources==0);
    SCITBX_ASSERT(num_fdblprime % diffBragg.sources == 0 );
    diffBragg.db_cryst.fpfdp_derivs.clear();
    //diffBragg.fdp_derivs.clear();
    int n_deriv = num_fprime / diffBragg.sources;
    for (int i_deriv=0; i_deriv < n_deriv; i_deriv++ ){
        for (int i=0; i<diffBragg.sources; i++){
            int idx=diffBragg.sources*i_deriv + i;
            double val = boost::python::extract<double>(fprime_derivs[idx]);
            diffBragg.db_cryst.fpfdp_derivs.push_back(val);
            val = boost::python::extract<double>(fdblprime_derivs[idx]);
            diffBragg.db_cryst.fpfdp_derivs.push_back(val);
        }
    }
  }


  static void set_fpfdp(simtbx::nanoBragg::diffBragg & diffBragg,
                boost::python::tuple const& fprime_fdblprime) {
    boost::python::list fprime = boost::python::extract<boost::python::list>(fprime_fdblprime[0]);
    boost::python::list fdblprime = boost::python::extract<boost::python::list>(fprime_fdblprime[1]);
    int num_fprime = boost::python::len(fprime);
    int num_fdblprime =   boost::python::len(fdblprime);
    SCITBX_ASSERT(num_fprime==diffBragg.sources);
    SCITBX_ASSERT(num_fdblprime==diffBragg.sources);
    diffBragg.db_cryst.fpfdp.clear();
    for (int i=0; i<diffBragg.sources; i++){
        double val = boost::python::extract<double>(fprime[i]);
        diffBragg.db_cryst.fpfdp.push_back(val);
        val = boost::python::extract<double>(fdblprime[i]);
        diffBragg.db_cryst.fpfdp.push_back(val);
    }
  }

  boost::python::tuple get_fpfdp(simtbx::nanoBragg::diffBragg const& diffBragg){
    // TODO implement
    boost::python::list fp;
    boost::python::list fdp;
    boost::python::tuple values;
    int n = diffBragg.db_cryst.fpfdp.size()/2;
    for(int i=0; i < n; i++ ){
        double fp_val = diffBragg.db_cryst.fpfdp[2*i];
        double fdp_val = diffBragg.db_cryst.fpfdp[2*i+1];
        fp.append(fp_val);
        fdp.append(fdp_val);
    }
    values = boost::python::make_tuple(fp,fdp);
    return values;
  }

  // change of basis operator see dxtbx.model.crystal.h
  static nanoBragg::mat3 get_O(simtbx::nanoBragg::diffBragg& diffBragg){
    return diffBragg.Omatrix;
  }

  void set_O(simtbx::nanoBragg::diffBragg& diffBragg, nanoBragg::mat3 const& value){
    diffBragg.Omatrix = value;
    diffBragg.db_cryst.eig_O << value[0], value[1], value[2],
            value[3], value[4], value[5],
            value[6], value[7], value[8];
    //std::cout << "eigO: " << diffBragg.eig_O << std::endl;
  }

  static nanoBragg::mat3 get_U(simtbx::nanoBragg::diffBragg& diffBragg){
    return diffBragg.Umatrix;
  }

  void set_U(simtbx::nanoBragg::diffBragg& diffBragg, nanoBragg::mat3 const& value){
    diffBragg.Umatrix = value;
    diffBragg.db_cryst.eig_U << value[0], value[1], value[2],
            value[3], value[4], value[5],
            value[6], value[7], value[8];
  }

  static nanoBragg::mat3 get_B(simtbx::nanoBragg::diffBragg& diffBragg){
    return diffBragg.Bmatrix;
  }

  void set_B(simtbx::nanoBragg::diffBragg& diffBragg, nanoBragg::mat3 const& value){
      diffBragg.Bmatrix = value;
      diffBragg.a_star[1] = value(0,0);  // initialize Angstrom version, update at the end
      diffBragg.a_star[2] = value(1,0);
      diffBragg.a_star[3] = value(2,0);
      diffBragg.b_star[1] = value(0,1);
      diffBragg.b_star[2] = value(1,1);
      diffBragg.b_star[3] = value(2,1);
      diffBragg.c_star[1] = value(0,2);
      diffBragg.c_star[2] = value(1,2);
      diffBragg.c_star[3] = value(2,2);

      diffBragg.user_cell = 0;
      diffBragg.user_matrix = 1;
      diffBragg.init_cell();

       diffBragg.db_cryst.eig_B << diffBragg.a0[1], diffBragg.b0[1], diffBragg.c0[1],
            diffBragg.a0[2], diffBragg.b0[2], diffBragg.c0[2],
            diffBragg.a0[3], diffBragg.b0[3], diffBragg.c0[3];
  }

  //void set_roi(simtbx::nanoBragg::diffBragg& diffBragg, boost::python::tuple const& tupleoftuples){
  //    boost::python::tuple frange;
  //    boost::python::tuple srange;
  //    frange = extract<boost::python::tuple>(tupleoftuples[0]);
  //    srange = extract<boost::python::tuple>(tupleoftuples[1]);
  //    diffBragg.roi_xmin=extract<int>(frange[0]);
  //    diffBragg.roi_xmax=extract<int>(frange[1]);
  //    diffBragg.roi_ymin=extract<int>(srange[0]);
  //    diffBragg.roi_ymax=extract<int>(srange[1]);
  //     /* fix this up so its not out-of-range */
  //    diffBragg.init_beamcenter();
  //    diffBragg.init_raw_pixels_roi();
  //    diffBragg.initialize_managers();
  //}

  /* region of interest in pixels */
  //static boost::python::tuple get_roi(simtbx::nanoBragg::diffBragg const& diffBragg) {
  //    boost::python::tuple frange;
  //    boost::python::tuple srange;
  //    frange = boost::python::make_tuple(diffBragg.roi_xmin,diffBragg.roi_xmax);
  //    srange = boost::python::make_tuple(diffBragg.roi_ymin,diffBragg.roi_ymax);
  //    return boost::python::make_tuple(frange,srange);
  //}

  static boost::python::tuple get_reference_origin(simtbx::nanoBragg::diffBragg & diffBragg){
    return boost::python::make_tuple(diffBragg.O_reference[0]*1000, diffBragg.O_reference[1]*1000, diffBragg.O_reference[2]*1000);
  }

  //simtbx::nanoBragg::af::shared<double> get_deriv_pix(simtbx::nanoBragg::diffBragg& diffBragg, int refine_id){
  //  simtbx::nanoBragg::af::shared<double> pixels = diffBragg.get_derivative_pixels(refine_id);
  //  // slice them
  //  //return pixels[simtbx::nanoBragg::af::slice(0,diffBragg.Npix_to_model)];
  //  return pixels[0:diffBragg.Npix_to_model];
  //}

  void set_reference_origin(simtbx::nanoBragg::diffBragg & diffBragg, boost::python::tuple const& values){
    diffBragg.O_reference[0] = extract<double>(values[0])/1000.;
    diffBragg.O_reference[1] = extract<double>(values[1])/1000.;
    diffBragg.O_reference[2] = extract<double>(values[2])/1000.;
  }

  void set_lambda_coef(simtbx::nanoBragg::diffBragg& diffBragg, boost::python::tuple const& values){
      double coef0 = extract<double>(values[0]);
      double coef1 = extract<double>(values[1]);
      diffBragg.lambda_managers[0]->value = coef0;
      diffBragg.lambda_managers[1]->value = coef1;
      diffBragg.use_lambda_coefficients = true;
  }

  static boost::python::tuple get_lambda_coef(simtbx::nanoBragg::diffBragg const& diffBragg) {
      double coef0=diffBragg.lambda_managers[0]->value;
      double coef1=diffBragg.lambda_managers[1]->value;
      return boost::python::make_tuple(coef0, coef1);
  }

  bool get_use_diffuse(simtbx::nanoBragg::diffBragg& diffBragg){
      return diffBragg.db_flags.use_diffuse;
  }

  void set_use_diffuse(simtbx::nanoBragg::diffBragg& diffBragg, bool val){
      diffBragg.db_flags.use_diffuse = val ;
  }

  static boost::python::tuple get_origin(simtbx::nanoBragg::diffBragg& diffBragg){
    return boost::python::make_tuple(diffBragg.pix0_vector[1]*1000,
                                    diffBragg.pix0_vector[2]*1000, diffBragg.pix0_vector[3]*1000);
  }

  static void  set_Fhkl_tuple(simtbx::nanoBragg::diffBragg& diffBragg, boost::python::tuple const& value) {
      //TODO nanoBragg set as well ?
      diffBragg.pythony_indices = extract<nanoBragg::indices >(value[0]);
      diffBragg.pythony_amplitudes = extract<nanoBragg::af::shared<double> >(value[1]);
      diffBragg.init_Fhkl();
      diffBragg.complex_miller = false;
      if (boost::python::len(value)==3){
          if (value[2]!=boost::python::object()){ // check if it is not None
              diffBragg.pythony_amplitudes2 = extract<nanoBragg::af::shared<double> >(value[2]);
              diffBragg.init_Fhkl2();
              diffBragg.complex_miller = true;
          }
      }
      diffBragg.linearize_Fhkl();
  }

  static boost::python::tuple get_Fhkl_tuple(simtbx::nanoBragg::diffBragg diffBragg) {
      int h,k,l;
      double temp;
      int hkls = diffBragg.h_range*diffBragg.k_range*diffBragg.l_range;
      diffBragg.pythony_indices = nanoBragg::indices(hkls,nanoBragg::af::init_functor_null<scitbx::vec3<int> >());
      diffBragg.pythony_amplitudes = nanoBragg::af::shared<double>(hkls,nanoBragg::af::init_functor_null<double>());
      diffBragg.pythony_amplitudes2 = nanoBragg::af::shared<double>(hkls,nanoBragg::af::init_functor_null<double>());
      int i=0;
      for(h=diffBragg.h_min;h<=diffBragg.h_max;++h){
          for(k=diffBragg.k_min;k<=diffBragg.k_max;++k){
              for(l=diffBragg.l_min;l<=diffBragg.l_max;++l){
                  if ( (h<=diffBragg.h_max) && (h>=diffBragg.h_min) && (k<=diffBragg.k_max)
                    && (k>=diffBragg.k_min) && (l<=diffBragg.l_max) && (l>=diffBragg.l_min)  ) {
                      /* populate all stored values */
                      nanoBragg::miller_t hkl (h,k,l);
                      diffBragg.pythony_indices[i] = hkl;
                      temp = diffBragg.Fhkl[h-diffBragg.h_min][k-diffBragg.k_min][l-diffBragg.l_min];
                      diffBragg.pythony_amplitudes[i] = temp;
                      if (diffBragg.complex_miller){
                          temp = diffBragg.Fhkl2[h-diffBragg.h_min][k-diffBragg.k_min][l-diffBragg.l_min];
                          diffBragg.pythony_amplitudes2[i] = temp;
                      }
                      ++i;
                  }
              }
          }
      }
      if (diffBragg.complex_miller)
        return boost::python::make_tuple(diffBragg.pythony_indices,diffBragg.pythony_amplitudes, diffBragg.pythony_amplitudes2);
      else
        return boost::python::make_tuple(diffBragg.pythony_indices,diffBragg.pythony_amplitudes);
  }

  void diffBragg_init_module() {
    using namespace boost::python;
    typedef return_value_policy<return_by_value> rbv;
    typedef default_call_policies dcp;
    typedef return_internal_reference<> rir;
    class_<simtbx::nanoBragg::diffBragg, bases<simtbx::nanoBragg::nanoBragg> >
            ("diffBragg", no_init)
      /* constructor that takes a dxtbx detector and beam model */
      .def(init<const dxtbx::model::Detector&,
                const dxtbx::model::Beam&,
                int >(
        (arg_("detector"),
         arg_("beam"),
         arg_("verbose")=0),
        "nanoBragg simulation initialized from dxtbx detector and beam objects"))

      /* function for doing some differentiation */
      //.def("add_diffBragg_spots",&simtbx::nanoBragg::diffBragg::add_diffBragg_spots,
      //  "gives derivative of average photon count w.r.t. parameter of choice")

      .def("add_diffBragg_spots", add_diffBragg_spots_A,
        "gives derivative of average photon count w.r.t. parameter of choice")

      .def("add_diffBragg_spots", add_diffBragg_spots_B,
        "gives derivative of average photon count w.r.t. parameter of choice")

      .def("add_diffBragg_spots", add_diffBragg_spots_C,
        "gives derivative of average photon count w.r.t. parameter of choice")

      /* Run this before running add_diffBragg_spots */
      .def("vectorize_umats",&simtbx::nanoBragg::diffBragg::vectorize_umats,
        "caches the UMATS")

      .def("initialize_managers",&simtbx::nanoBragg::diffBragg::initialize_managers,
        "sets up the managers")

      .def("refine", &simtbx::nanoBragg::diffBragg::refine,
        "sets the refinement flag to true")

      .def("fix", &simtbx::nanoBragg::diffBragg::fix,
        "sets the refinement flag to false")

      .def("let_loose", &simtbx::nanoBragg::diffBragg::let_loose,
        "sets the refinement flag to True(but doesnt initialize)")

      .def("get_value", &simtbx::nanoBragg::diffBragg::get_value, "get value of the refinement parameter")

      .def("set_value", &simtbx::nanoBragg::diffBragg::set_value, "set value of the refinement parameter")

      .def("set_ncells_values", &simtbx::nanoBragg::diffBragg::set_ncells_values, "set Ncells values as a 3-tuple (Na, Nb, Nc)")

      .def("get_ncells_values", &simtbx::nanoBragg::diffBragg::get_ncells_values, "get Ncells values as a 3-tuple (Na, Nb, Nc)")

      .def("add_diffBragg_spots_full", &simtbx::nanoBragg::diffBragg::add_diffBragg_spots_full, "forward model and gradients at every pixel")

      .def("get_ncells_derivative_pixels", &simtbx::nanoBragg::diffBragg::get_ncells_derivative_pixels, "get derivatives of intensity w.r.t (Na, Nb, Nc)")

      .def("get_fp_fdp_derivative_pixels", &simtbx::nanoBragg::diffBragg::get_fp_fdp_derivative_pixels, "get derivatives of intensity w.r.t c,d that describe fprime and fdblprime (see diffBragg.utils)")

      .def("get_ncells_def_derivative_pixels", &simtbx::nanoBragg::diffBragg::get_ncells_def_derivative_pixels, "get derivatives of intensity w.r.t (Nd, Ne, Nf)")

      //.def("get_Na_derivative_pixels", &simtbx::nanoBragg::diffBragg::get_Na_derivative_pixels, "get derivatives of intensity w.r.t (Nd, Ne, Nf)")

      .def("get_aniso_eta_deriv_pixels", &simtbx::nanoBragg::diffBragg::get_aniso_eta_deriv_pixels, "get derivatives of intensity w.r.t anisotropic mosaicity model")

      .def("get_aniso_eta_second_deriv_pixels", &simtbx::nanoBragg::diffBragg::get_aniso_eta_second_deriv_pixels, "get second derivatives of intensity w.r.t anisotropic mosaicity model")

      .def("get_ncells_def_second_derivative_pixels", &simtbx::nanoBragg::diffBragg::get_ncells_def_second_derivative_pixels, "get 2nd derivatives of intensity w.r.t (Nd, Ne, Nf)")

      .def("get_ncells_second_derivative_pixels", &simtbx::nanoBragg::diffBragg::get_ncells_second_derivative_pixels, "get second derivatives of intensity w.r.t (Na, Nb, Nc)")

      .def("get_lambda_derivative_pixels", &simtbx::nanoBragg::diffBragg::get_lambda_derivative_pixels, "get derivatives of intensity w.r.t (lambda0, lambda1) as a 2-tuple")

      .def("print_if_refining", &simtbx::nanoBragg::diffBragg::print_if_refining, "print a statement if the refiner flag is set")

      .def("quick_Fhkl_update", &simtbx::nanoBragg::diffBragg::quick_Fcell_update, "update Fcell amplitudes without reallocating memory")
      .def("update_xray_beams", &simtbx::nanoBragg::diffBragg::update_xray_beams, "update xray beams without reallocating memory")

      .def("set_ucell_derivative_matrix",  &simtbx::nanoBragg::diffBragg::set_ucell_derivative_matrix, "Boo-ya")

      .def("set_ucell_second_derivative_matrix",  &simtbx::nanoBragg::diffBragg::set_ucell_second_derivative_matrix,
            "Two-ya")

      //.def("get_ucell_derivative_matrix",  &simtbx::nanoBragg::diffBragg::get_ucell_derivative_matrix, "scooby snacks")

      //.def("get_derivative_pixels", get_deriv_pix,
      //      "gets the manager raw image containing first derivatives")

      .def("__get_derivative_pixels", &simtbx::nanoBragg::diffBragg::get_derivative_pixels,
            "gets the manager raw image containing first derivatives")

      .def("get_second_derivative_pixels", &simtbx::nanoBragg::diffBragg::get_second_derivative_pixels,
            "gets the manager raw image containing second derivatives")

      .def("zero_raw_pixel_rois", &simtbx::nanoBragg::diffBragg::zero_raw_pixel_rois,
           "reallocate the raw image ROIs (they are usually tiny)")

      .def("get_origin", get_origin, "get panel origin")

      .def("shift_origin_z",
            &simtbx::nanoBragg::diffBragg::shift_originZ,
            "shift the Z-component of the pixel0 vector in every detector panel")

      .def("update_dxtbx_geoms",
            &simtbx::nanoBragg::diffBragg::update_dxtbx_geoms,
             (arg_("detector"),
              arg_("beam"),
              arg_("panel_id")=0,
              arg_("panel_rot_angO")=0,
              arg_("panel_rot_angF")=0,
              arg_("panel_rot_angS")=0,
              arg_("panel_offsetX")=0,
              arg_("panel_offsetY")=0,
              arg_("panel_offsetZ")=0, arg_("force")=true),
           "update the geometries with new dxtbx models, number of pixels should remain constant")

      .def("free_Fhkl2",&simtbx::nanoBragg::diffBragg::free_Fhkl2)

      .def("show_fp_fdp", &simtbx::nanoBragg::diffBragg::show_fp_fdp)

      .def("show_heavy_atom_data", &simtbx::nanoBragg::diffBragg::show_heavy_atom_data)

#ifdef NANOBRAGG_HAVE_CUDA
      .def("gpu_free",&simtbx::nanoBragg::diffBragg::gpu_free)
#endif

      .def("set_mosaic_blocks_prime",
             &simtbx::nanoBragg::diffBragg::set_mosaic_blocks_prime,
             "enter a list of unitary matrix derivatives Uprime")

      .def("set_mosaic_blocks_dbl_prime",
             &simtbx::nanoBragg::diffBragg::set_mosaic_blocks_dbl_prime,
             "enter a list of unitary matrix 2nd derivatives Uprime")

      .def("show_timings",
           &simtbx::nanoBragg::diffBragg::show_timing_stats,
           (arg_("MPI_RANK")=0),
           "Shows a few key accumulated timings (in milliseconds) of specific code portions")

      //.add_property("region_of_interest",
      //       make_function(&get_roi,rbv()),
      //       make_function(&set_roi,dcp()),
      //       "region of interest on detector: fast_min fast_max slow_min slow_max")

      .add_property("raw_pixels_roi",
                     make_getter(&simtbx::nanoBragg::diffBragg::raw_pixels_roi,rbv()),
                     make_setter(&simtbx::nanoBragg::diffBragg::raw_pixels_roi,dcp()),
                    "raw pixels from region of interest only")

      .add_property("oversample_omega",
                     make_getter(&simtbx::nanoBragg::diffBragg::oversample_omega,rbv()),
                     make_setter(&simtbx::nanoBragg::diffBragg::oversample_omega,dcp()),
                    "whether to use an average solid angle correction per pixel, or one at the sub pixel level")

      .add_property("track_Fhkl",
                     make_getter(&simtbx::nanoBragg::diffBragg::track_Fhkl,rbv()),
                     make_setter(&simtbx::nanoBragg::diffBragg::track_Fhkl,dcp()),
                    "whether to track Fhkl per pixel")

      .add_property("fprime_fdblprime",
             make_function(&get_fpfdp,rbv()),
             make_function(&set_fpfdp,dcp()),
             "fprime and fdoubleprime, passed as two lists whose length should be the same")

      .add_property("fprime_fdblprime_derivs",
             make_function(&get_fpfdp_derivs,rbv()), // TODO implement me
             make_function(&set_fpfdp_derivs,dcp()),
             "fprime and fdoubleprime derivatives, as a 2-tuple. e.g. 1st tuple element is fprime derivs, length should be multiple of number of sources")

      .add_property("heavy_atom_data",
             make_function(&get_atom_data,rbv()),
             make_function(&set_atom_data,dcp()),
             "5 tuple of fracx_list,fracy_list,fracz_list,Bfactor_list,Occupancy_list, each list is same length")

      .add_property("no_Nabc_scale",
                     make_getter(&simtbx::nanoBragg::diffBragg::no_Nabc_scale,rbv()),
                     make_setter(&simtbx::nanoBragg::diffBragg::no_Nabc_scale,dcp()),
                    "toggle off the Nabc scale factor in the forward model (such that it is replaced entirely by spot_scale)")

      .add_property("__number_of_pixels_modeled_using_diffBragg", // protect this by making it a long name
                     make_getter(&simtbx::nanoBragg::diffBragg::Npix_to_model,rbv()),
                     make_setter(&simtbx::nanoBragg::diffBragg::Npix_to_model,dcp()),
                    "num pix modeled most recently")

      .add_property("compute_curvatures",
                     make_getter(&simtbx::nanoBragg::diffBragg::compute_curvatures,rbv()),
                     make_setter(&simtbx::nanoBragg::diffBragg::compute_curvatures,dcp()),
                    "Whether to compute the curvatures")

      .add_property("only_save_omega_kahn",
                     make_getter(&simtbx::nanoBragg::diffBragg::only_save_omega_kahn,rbv()),
                     make_setter(&simtbx::nanoBragg::diffBragg::only_save_omega_kahn,dcp()),
                    "ONLY simulate the kahn polarization correction and solid angle components")

      .add_property("reference_origin",
            make_function(&get_reference_origin, rbv()),
            make_function(&set_reference_origin, dcp()),
            "reference origin for doing panel rotations for a group of panels")

      .add_property("Umatrix",
             make_function(&get_U,rbv()),
             make_function(&set_U,dcp()),
             "rotation portion of Amatrix (dxtbx format, C.get_U())")

      .add_property("Bmatrix",
             make_function(&get_B,rbv()),
             make_function(&set_B,dcp()),
             "unit cell portion of Amatrix (dxtbx format; C.get_B())")

      .add_property("Omatrix",
             make_function(&get_O,rbv()),
             make_function(&set_O,dcp()),
             "change of basis Omatrix")

      .add_property("Ncells_abc",
             make_function(&get_Nabc,rbv()),
             make_function(&set_Nabc,dcp()),
             "number of mosaic domains along an axis")

      .add_property("Ncells_def",
             make_function(&get_Ndef,rbv()),
             make_function(&set_Ndef,dcp()),
             "off-axis components in mosaic domain size matrix")

      .add_property("Ncells_abc_aniso",
             make_function(&get_Nabc_aniso,rbv()),
             make_function(&set_Nabc_aniso,dcp()),
             "number of mosaic domains along an axis, can be anisotropic")

      .add_property("update_oversample_during_refinement",
                     make_getter(&simtbx::nanoBragg::diffBragg::update_oversample_during_refinement,rbv()),
                     make_setter(&simtbx::nanoBragg::diffBragg::update_oversample_during_refinement,dcp()),
                     "Allows oversample to change as Ncells abc changes")

      .add_property("isotropic_ncells",
                     make_getter(&simtbx::nanoBragg::diffBragg::isotropic_ncells,rbv()),
                     make_setter(&simtbx::nanoBragg::diffBragg::isotropic_ncells,dcp()),
                     "refine (Na, Nb, Nc) as three independent parameters")

      .add_property("has_anisotropic_mosaic_spread",
                     make_getter(&simtbx::nanoBragg::diffBragg::modeling_anisotropic_mosaic_spread,rbv()),
                     make_setter(&simtbx::nanoBragg::diffBragg::modeling_anisotropic_mosaic_spread,dcp()),
                     "compute derivatives for a single mosaic spread parameter")

      .add_property("use_lambda_coefficients",
                     make_getter(&simtbx::nanoBragg::diffBragg::use_lambda_coefficients,rbv()),
                     make_setter(&simtbx::nanoBragg::diffBragg::use_lambda_coefficients,dcp()),
                     "represent wavelength using coefficients a + b*source where a,b are coefficients and source is the source index")

      .add_property("Fhkl_tuple",
            make_function(&get_Fhkl_tuple,rbv()),
            make_function(&set_Fhkl_tuple,dcp()),
            "hkl and Freal and Fcomplex as a 3-tuple of (indices ,flex-double, flex-double). experimental, also accepts (indices, flex-double, None) for mod only")

      .add_property("lambda_coefficients",
             make_function(&get_lambda_coef,rbv()),
             make_function(&set_lambda_coef,dcp()),
             "coefficients for source_lambda refinement: `lambda = coef0 + coef1*source`  where `source` is the source index")

     // CUDA PROPERTIES
      .add_property("use_cuda",
             make_getter(&simtbx::nanoBragg::diffBragg::use_cuda,rbv()),
             make_setter(&simtbx::nanoBragg::diffBragg::use_cuda,dcp()),
             "use GPU acceleration")

      .add_property("record_time",
             make_getter(&simtbx::nanoBragg::diffBragg::record_time,rbv()),
             make_setter(&simtbx::nanoBragg::diffBragg::record_time,dcp()),
             "store timing data for key aspects of the code")

      .add_property("Npix_to_allocate",
            make_getter(&simtbx::nanoBragg::diffBragg::Npix_to_allocate,rbv()),
            make_setter(&simtbx::nanoBragg::diffBragg::Npix_to_allocate,dcp()),
            "number of pixels per parameter to allocate in the GPU (default value of -1 means auto-determine)")

      .add_property("use_diffuse",
            make_function(&get_use_diffuse,rbv()),
            make_function(&set_use_diffuse,dcp()),
            "sim with diffuse")

    ; // end of diffBragg extention

  } // end of diffBragg_init_module
} // end of namespace
} // end of boost python namespace
} // end of refine
} // end of simtbx

BOOST_PYTHON_MODULE(simtbx_diffBragg_ext)
{
  simtbx::diffBragg::boost_python::diffBragg_init_module();
}
