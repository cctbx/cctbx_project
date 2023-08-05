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

#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <scitbx/array_family/boost_python/selections_wrapper.h>

using namespace boost::python;
using namespace molprobity::reduce;

BOOST_PYTHON_MODULE(mmtbx_reduce_ext)
{
  // Dependencies

  // Describe and name compound classes that we need access to beyond those that are
  // already defined for us by scitbx arrays that are defined elsewhere.
  typedef scitbx::af::boost_python::shared_wrapper<double> wdbl;
  class_<wdbl::w_t> wd = wdbl::wrap("af_shared_double");
  scitbx::af::boost_python::select_wrappers<
    double, scitbx::af::shared<double> >::wrap(wd);

  typedef scitbx::af::shared<bool> afsbool;
  typedef scitbx::af::boost_python::shared_wrapper<afsbool> wwbool;
  class_<wwbool::w_t> wwb = wwbool::wrap("af_shared_af_shared_bool");
  scitbx::af::boost_python::select_wrappers<
    afsbool, scitbx::af::shared<afsbool> >::wrap(wwb);

  typedef scitbx::af::shared<molprobity::probe::Point> afsPoint;
  typedef scitbx::af::boost_python::shared_wrapper<afsPoint> wwPoint;
  class_<wwPoint::w_t> wwd = wwPoint::wrap("af_shared_af_shared_Point");
  scitbx::af::boost_python::select_wrappers<
    afsPoint, scitbx::af::shared<afsPoint> >::wrap(wwd);

  typedef scitbx::af::shared<molprobity::probe::ExtraAtomInfo> afsei;
  typedef scitbx::af::boost_python::shared_wrapper<afsei> wwei;
  class_<wwei::w_t> wwExtraInfo = wwei::wrap("af_shared_af_shared_ExtraAtomInfo");
  scitbx::af::boost_python::select_wrappers<
    afsei, scitbx::af::shared<afsei> >::wrap(wwExtraInfo);

  // Define the flex array wrapping for these classes because we take them as parameters.
  scitbx::boost_python::container_conversions::tuple_mapping_variable_capacity<
    scitbx::af::shared< scitbx::af::shared<molprobity::probe::Point> > >();
  scitbx::boost_python::container_conversions::tuple_mapping_variable_capacity<
    scitbx::af::shared< scitbx::af::shared<molprobity::probe::ExtraAtomInfo> > >();
  scitbx::boost_python::container_conversions::tuple_mapping_variable_capacity<
    scitbx::af::shared< scitbx::af::shared<bool> > >();

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
}
