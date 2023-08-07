// Copyright(c) 2023, Richardson Lab at Duke
// Licensed under the Apache 2 license
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissionsand
// limitations under the License.

// Enable functions with up to 20 parameters to be called.  Default of 15 is insufficient
#define BOOST_PYTHON_MAX_ARITY 20
#include <boost/python.hpp>
//#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <scitbx/boost_python/container_conversions.h>
#include <mmtbx/reduce/PositionReturn.h>
#include <mmtbx/reduce/InteractionGraph.h>

#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/selections_wrapper.h>

using namespace boost::python;
using namespace molprobity::reduce;

BOOST_PYTHON_MODULE(mmtbx_reduce_ext)
{
  // Dependencies
  //boost::python::object pdb_hierarchy_ext = boost::python::import( "iotbx_pdb_hierarchy_ext" );
  //boost::python::object pdb_ext = boost::python::import( "iotbx_pdb_ext" );
  //boost::python::object rstbx_ext = boost::python::import("rstbx_array_family_flex_ext");
  //boost::python::object scitbx_flex_ext = boost::python::import("scitbx_array_family_flex_ext");
  //boost::python::object scitbx_shared_ext = boost::python::import("scitbx_array_family_shared_ext");

  // Describe and name compound classes that we need access to beyond those that are
  // already defined for us by scitbx arrays that are defined elsewhere.

  /*
  std::cout << "XXX Mapping the double" << std::endl;
  typedef scitbx::af::boost_python::shared_wrapper<double> wdbl;
  class_<wdbl::w_t> wd = wdbl::wrap("af_shared_double");
  scitbx::af::boost_python::select_wrappers<
    double, scitbx::af::shared<double> >::wrap(wd);
  std::cout << "XXX Done mapping the double" << std::endl;
  /// @todo This does not make the type available to Python
  //scitbx::af::boost_python::flex_wrapper<scitbx::af::shared<double>>::plain("flex_double");
  */

  typedef scitbx::af::shared<bool> afsbool;
  typedef scitbx::af::boost_python::shared_wrapper<afsbool> wwbool;
  class_<wwbool::w_t> wwb = wwbool::wrap("af_shared_af_shared_bool");
  scitbx::af::boost_python::select_wrappers<
    afsbool, scitbx::af::shared<afsbool> >::wrap(wwb);

  /*
  std::cout << "XXX Mapping the Point" << std::endl;
  typedef scitbx::af::shared<molprobity::probe::Point> afsPoint;
  typedef scitbx::af::boost_python::shared_wrapper<afsPoint> wwPoint;
  class_<wwPoint::w_t> wwd = wwPoint::wrap("af_shared_af_shared_Point");
  scitbx::af::boost_python::select_wrappers<
    afsPoint, scitbx::af::shared<afsPoint> >::wrap(wwd);
  std::cout << "XXX Done mapping the Point" << std::endl;
  */

  typedef scitbx::af::shared<molprobity::probe::ExtraAtomInfo> afsei;
  typedef scitbx::af::boost_python::shared_wrapper<afsei> wwei;
  class_<wwei::w_t> wwExtraInfo = wwei::wrap("af_shared_af_shared_ExtraAtomInfo");
  scitbx::af::boost_python::select_wrappers<
    afsei, scitbx::af::shared<afsei> >::wrap(wwExtraInfo);

  // Define the flex array wrapping for these classes because we take them as parameters.
  scitbx::boost_python::container_conversions::tuple_mapping_variable_capacity<
    scitbx::af::shared< scitbx::af::shared<molprobity::probe::Point> > >();
  /*
  scitbx::boost_python::container_conversions::tuple_mapping_variable_capacity<
    scitbx::af::shared< scitbx::af::shared<molprobity::probe::ExtraAtomInfo> > >();
  scitbx::boost_python::container_conversions::tuple_mapping_variable_capacity<
    scitbx::af::shared< scitbx::af::shared<bool> > >();
    */

  class_<PositionReturn>("PositionReturn")
    .def(init<>())
    .def(init< scitbx::af::shared<iotbx::pdb::hierarchy::atom>
      , scitbx::af::shared< scitbx::af::shared<molprobity::probe::Point> >
      , scitbx::af::shared< scitbx::af::shared<molprobity::probe::ExtraAtomInfo> >
      , scitbx::af::shared< scitbx::af::shared<bool> >
      , scitbx::af::shared<double>
      >()
    )
    .add_property("atoms", &PositionReturn::atoms)
    .add_property("positions", &PositionReturn::positions)
    .add_property("extraInfos", &PositionReturn::extraInfos)
    .add_property("deleteMes", &PositionReturn::deleteMes)
    .add_property("preferenceEnergies", &PositionReturn::preferenceEnergies)
    ;
  // Export the global functions
  def("PositionReturn_test", PositionReturn_test, "Test all classes defined in PositionReturn.h.");

  def("PairsOverlap", PairsOverlap, "Test for overlap between two pairs of atoms.");
}
