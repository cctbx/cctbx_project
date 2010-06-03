#include <boost/python.hpp>
#include <spotfinder/core_toolbox/distl.h>
#include <spotfinder/core_toolbox/spotfilter.h>
#include <spotfinder/core_toolbox/boost_python/singlemask.h>
#include <scitbx/array_family/flex_types.h>

using namespace boost::python;
using namespace spotfinder::distltbx;
namespace af = scitbx::af;

namespace{

af::shared<int>
bin_populations(af::const_ref<double> const& resolution, double const& limit){
  int Total_rows = int(std::pow(limit/resolution[resolution.size()-1],3.))+1;
  af::shared<double> BinRes;
  for (int i=0; i<Total_rows; ++i) {
    BinRes.push_back(limit/std::pow(i+1.,1.0/3.0));
  }
  af::shared<int> z(Total_rows);
  int row=0;
  for (const double* resptr = resolution.begin();
       resptr != resolution.end();++resptr){
    if (*resptr >= BinRes[row]) {
      z[row]+=1;
    } else {
      while (*resptr < BinRes[row+1]) {
        row+=1;
      }
      row+=1;
      z[row]=1;
    }
  }
  return z;
}

af::shared<int>
bin_exclusion_and_impact(af::flex_double const& sorted_resolutions,
                       af::flex_int const& sorted_order,
                       af::flex_double const& limits_low_hi_res,
                       af::flex_int & bin_impact_output_array){
    af::shared<int>sorted_order_subset;
    for (int i = 0; i<sorted_resolutions.size(); ++i){
      bool clear = true;
      for (int s = 0; s < limits_low_hi_res.size(); s+=2){
        if ( (limits_low_hi_res[s+1] < sorted_resolutions[i]) &&
             (sorted_resolutions[i] < limits_low_hi_res[s]) ){
          clear=false;
          bin_impact_output_array[s/2]+=1;
          break;
        }
      }
      if (clear){ sorted_order_subset.push_back(sorted_order[i]); }
    }
    return sorted_order_subset;
}

}

namespace spotfinder { namespace distltbx { namespace boost_python {

  void wrap_geometry_2d();

namespace {

  void init_module() {
    using namespace boost::python;
    wrap_geometry_2d();
  }

} // namespace <anonymous>
}}} // namespace spotfinder::distltbx::boost_python


BOOST_PYTHON_MODULE(spotfinder_distltbx_ext)
{
   typedef return_value_policy<return_by_value> rbv;
   typedef default_call_policies dcp;

   spotfinder::distltbx::boost_python::init_module();

   class_<Distl::point>("distl_point",init<int, int>())
     .def_readonly("x",             &Distl::point::x)
     .def_readonly("y",             &Distl::point::y)
   ;

   class_<Distl::icering>("distl_icering", no_init)
     .add_property("lowerr2",make_getter(&Distl::icering::lowerr2,rbv()))
     .add_property("upperr2",make_getter(&Distl::icering::upperr2,rbv()))
     .add_property("lowerresol",make_getter(&Distl::icering::lowerresol,rbv()))
     .add_property("upperresol",make_getter(&Distl::icering::upperresol,rbv()))
     .add_property("strength",make_getter(&Distl::icering::strength,rbv()))
   ;

   class_<Distl::spot_base >("spot_base", no_init)
     .def("area",                  &Distl::spot_base::area)
     .add_property("bodypixels",make_getter(&Distl::spot_base::bodypixels,rbv()),
                                make_setter(&Distl::spot_base::bodypixels,dcp()))
     .add_property("peak",make_getter(&Distl::spot_base::peak,rbv()),
                          make_setter(&Distl::spot_base::peak,dcp()))
     .add_property("maximas",make_getter(&Distl::spot_base::maximas,rbv()))
   ;
   class_<Distl::spot_shapes >("spot_shapes", no_init)
     .def("ctr_mass_x",                  &Distl::spot_shapes::ctr_mass_x)
     .def("ctr_mass_y",                  &Distl::spot_shapes::ctr_mass_y)
     .def("eigenvalue",                  &Distl::spot_shapes::eigenvalue)
     .def("eigenvector",                 &Distl::spot_shapes::eigenvector)
     .def("model_eccentricity",  &Distl::spot_shapes::model_eccentricity)
     .def("a",                           &Distl::spot_shapes::a)
     .def("b",                           &Distl::spot_shapes::b)
     .def("model_ellipse",               &Distl::spot_shapes::model_ellipse)
     .def("show_axes",                   &Distl::spot_shapes::show_axes)
     .def("com_valid",                   &Distl::spot_shapes::com_valid)
     .add_property("total_mass",make_getter(&Distl::spot_shapes::total_mass,rbv()))
   ;
   class_<w_spot, bases<Distl::spot_base,Distl::spot_shapes > >
     ("distl_spot", init<>())
     .def("max_pxl_x",                     &w_spot::max_pxl_x)
     .def("max_pxl_y",                     &w_spot::max_pxl_y)
     .def_readonly("nmaxima",      &w_spot::nmaxima)
     .def("intensity",             &w_spot::intensity)
     .def("perimeter",             &w_spot::perimeter)
     .def("shape",                 &w_spot::get_shape)
     .def("majoraxis",             &w_spot::get_majoraxis)
     .def("minoraxis",             &w_spot::get_minoraxis)
     .def("skewness",              &w_spot::skewness)
     .add_property("peakintensity",make_getter(&w_spot::peakintensity,rbv()),
                                   make_setter(&w_spot::peakintensity,dcp()))
     .add_property("wts",          make_getter(&w_spot::wts,rbv()),
                                   make_setter(&w_spot::wts,dcp()))
     //type scitbx::sf::shared must be exposed using add_property,
     // which can specify a return value policy, rather than def_readonly
   ;

   class_<w_Distl>("w_Distl", init<std::string,bool>())
     .def("set_resolution_outer",&w_Distl::set_resolution_outer)
     .def("setspotimg",&w_Distl::setspotimg)
     .def("set_tiling",&w_Distl::set_tiling)
     .def("Z_data",&w_Distl::Z_data)
     .def("mod_data",&w_Distl::mod_data)
     .def("get_underload",&w_Distl::get_underload)
     .def("pxlclassify",&w_Distl::pxlclassify)
     .def("search_icerings",&w_Distl::search_icerings)
     .def("search_maximas",&w_Distl::search_maximas)
     .def("search_spots",&w_Distl::search_spots)
     .def("search_overloadpatches",&w_Distl::search_overloadpatches)
     .def("finish_analysis",&w_Distl::finish_analysis)
     .add_property("spots",make_getter(&w_Distl::spots,rbv()))

     .def("nicerings",&w_Distl::nicerings)
     .def("isIsolated",&w_Distl::isIsolated)
     .add_property("icerings",make_getter(&w_Distl::icerings,rbv()))
     .def("imgresol",&w_Distl::imgresol)
   ;

   def("linear_char_scratchpad",&spotfinder::distltbx::linear_char_scratchpad);
   def("bin_populations",bin_populations);
   def("bin_exclusion_and_impact",bin_exclusion_and_impact);

   using boost::python::arg;
   class_<SpotFilterAgent>("SpotFilterAgent", init< double const&,
      double const&, double const&, double const&, double const&,
      af::shared<Distl::icering> >(
        (arg("pixel_size"),arg("xbeam"),arg("ybeam"),arg("distance"),
         arg("wavelength"),arg("icerings")
        )
      ))
      .def("precompute_resolution",&SpotFilterAgent::precompute_resolution)
      .def("filter",&SpotFilterAgent::filter)
      .def("resolution_sort",&SpotFilterAgent::resolution_sort)
      .def("resolution_sort_nztt",&SpotFilterAgent::resolution_sort_nztt)
      .def("get_resolution",&SpotFilterAgent::get_resolution)
      .def("get_property",&SpotFilterAgent::get_property)
      .def("set_arguments",&SpotFilterAgent::set_arguments)
      .def("order_by",&SpotFilterAgent::order_by)
      .enable_pickling()
      .def_readonly("pixel_size",&SpotFilterAgent::pixel_size)
      .def_readonly("xbeam",&SpotFilterAgent::xbeam)
      .def_readonly("ybeam",&SpotFilterAgent::ybeam)
      .def_readonly("distance",&SpotFilterAgent::distance)
      .def_readonly("wavelength",&SpotFilterAgent::wavelength)
      .add_property("icerings",make_getter(&SpotFilterAgent::icerings,rbv()))
   ;
   class_<SingleMask>("SingleMask", init< af::shared<Distl::spot>,
      af::shared<int> >()
      )
      .def_readonly("x",&SingleMask::x)
      .def_readonly("y",&SingleMask::y)
   ;
}
