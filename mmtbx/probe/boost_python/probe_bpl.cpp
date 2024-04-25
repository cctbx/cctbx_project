// Copyright(c) 2021, Richardson Lab at Duke
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
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <scitbx/boost_python/container_conversions.h>
#include <mmtbx/probe/DotSpheres.h>
#include <mmtbx/probe/SpatialQuery.h>
#include <mmtbx/probe/Scoring.h>

using namespace boost::python;
using namespace molprobity::probe;

/// @brief Helper function to wrap the Point internal array so we can read its elements.
///
/// If you want to change values in a Point, construct a new one using the 3-parameter
/// constructor.
boost::python::tuple wrap_vec3_array(Point const& d) {
  boost::python::list a;
  for (size_t i = 0; i < d.size(); ++i) {
    a.append(d.elems[i]);
  }
  return boost::python::tuple(a);
}

/// @brief Helper function to set the i_seq on an atom.
void set_atom_i_seq(iotbx::pdb::hierarchy::atom& atom, int i_seq) {
  atom.data->i_seq = i_seq;
}

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(check_dot_overloads, DotScorer::check_dot, 5, 6)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(interaction_type_overloads, DotScorer::interaction_type, 2, 3)

BOOST_PYTHON_MODULE(mmtbx_probe_ext)
{
  // Dependencies
  //boost::python::object pdb_hierarchy_ext = boost::python::import( "iotbx_pdb_hierarchy_ext" );
  //boost::python::object pdb_ext = boost::python::import( "iotbx_pdb_ext" );

  // Describe and name compound classes that we need access to beyond those that are
  // already defined for us by scitbx arrays that are defined elsewhere.

  class_<ContactResult>("ContactResult", init<>())
    .add_property("closestContact", &ContactResult::closestContact)
    .add_property("distAboveSurface", &ContactResult::distAboveSurface)
    ;

  class_<ExtraAtomInfo>("ExtraAtomInfo")
    .def(init< optional<double, bool, bool, bool, bool, int, std::string> >())
    .def(init<ExtraAtomInfo const &>())
    .add_property("vdwRadius", &ExtraAtomInfo::getVdwRadius, &ExtraAtomInfo::setVdwRadius)
    .add_property("isAcceptor", &ExtraAtomInfo::getIsAcceptor, &ExtraAtomInfo::setIsAcceptor)
    .add_property("isDonor", &ExtraAtomInfo::getIsDonor, &ExtraAtomInfo::setIsDonor)
    .add_property("isDummyHydrogen", &ExtraAtomInfo::getIsDummyHydrogen, &ExtraAtomInfo::setIsDummyHydrogen)
    .add_property("isIon", &ExtraAtomInfo::getIsIon, &ExtraAtomInfo::setIsIon)
    .add_property("charge", &ExtraAtomInfo::getCharge, &ExtraAtomInfo::setCharge)
    .add_property("altLoc", &ExtraAtomInfo::getAltLoc, &ExtraAtomInfo::setAltLoc)
    ;
  // Define the flex array wrapping for this class because we take it as a parameter.
  scitbx::boost_python::container_conversions::tuple_mapping_variable_capacity<
    scitbx::af::shared<ExtraAtomInfo> >();

  class_<ExtraAtomInfoMap>("ExtraAtomInfoMap", init< scitbx::af::shared<iotbx::pdb::hierarchy::atom>, scitbx::af::shared<ExtraAtomInfo> >())
    .def("getMappingFor", &ExtraAtomInfoMap::getMappingFor, return_internal_reference<>())
    .def("setMappingFor", &ExtraAtomInfoMap::setMappingFor)
    ;

  enum_<DotScorer::OverlapType>("OverlapType")
    .value("Ignore", DotScorer::Ignore)
    .value("NoOverlap", DotScorer::NoOverlap)
    .value("Clash", DotScorer::Clash)
    .value("HydrogenBond", DotScorer::HydrogenBond)
    ;

  enum_<DotScorer::InteractionType>("InteractionType")
    .value("WideContact", DotScorer::WideContact)
    .value("CloseContact", DotScorer::CloseContact)
    .value("WeakHydrogenBond", DotScorer::WeakHydrogenBond)
    .value("SmallOverlap", DotScorer::SmallOverlap)
    .value("Bump", DotScorer::Bump)
    .value("BadBump", DotScorer::BadBump)
    .value("StandardHydrogenBond", DotScorer::StandardHydrogenBond)
    .value("Invalid", DotScorer::Invalid)
    ;

  class_<DotScorer::CheckDotResult>("CheckDotResult", init<>())
    .add_property("overlapType", &DotScorer::CheckDotResult::overlapType)
    .add_property("cause", &DotScorer::CheckDotResult::cause)
    .add_property("overlap", &DotScorer::CheckDotResult::overlap)
    .add_property("gap", &DotScorer::CheckDotResult::gap)
    .add_property("annular", &DotScorer::CheckDotResult::annular)
    ;

  class_<DotScorer::ScoreDotsResult>("ScoreDotsResult", init<>())
    .add_property("valid", &DotScorer::ScoreDotsResult::valid)
    .add_property("bumpSubScore", &DotScorer::ScoreDotsResult::bumpSubScore)
    .add_property("hBondSubScore", &DotScorer::ScoreDotsResult::hBondSubScore)
    .add_property("attractSubScore", &DotScorer::ScoreDotsResult::attractSubScore)
    .add_property("hasBadBump", &DotScorer::ScoreDotsResult::hasBadBump)
    .def("totalScore", &DotScorer::ScoreDotsResult::totalScore)
    ;

  class_<DotSphere>("DotSphere", init<double, double>())
    .def(init<>())
    .def("dots", &DotSphere::dotsCopyForPythonWrapping)
    .def("radius", &DotSphere::radius)
    .def("density", &DotSphere::density)
    .def("test", &DotSphere::test)
  ;

  class_<DotSphereCache>("DotSphereCache", init<double>())
    .def("get_sphere", &DotSphereCache::get_sphere, return_internal_reference<>())
    .def("size", &DotSphereCache::size)
    .def("test", &DotSphereCache::test)
  ;

  class_<SpatialQuery>("SpatialQuery", init<Point, Point, Point>())
    .def(init<scitbx::af::shared<iotbx::pdb::hierarchy::atom> const>())
    .def("add", &SpatialQuery::add)
    .def("remove", &SpatialQuery::remove)
    .def("neighbors", &SpatialQuery::neighbors)
    .def("test", &SpatialQuery::test)
  ;

  class_<DotScorer>("DotScorer",
        init< ExtraAtomInfoMap,
        optional<double, double, double, double, double, double, double, double,bool,bool> >())
    .def("point_inside_atoms", &DotScorer::point_inside_atoms)
    .def("trim_dots", &DotScorer::trim_dots)
    .def("check_dot", &DotScorer::check_dot, check_dot_overloads())
    .def("count_surface_dots", &DotScorer::count_surface_dots)
    .def("score_dots", &DotScorer::score_dots)
    .def("interaction_type", &DotScorer::interaction_type, interaction_type_overloads())
    .def("interaction_type_name", &DotScorer::interaction_type_name)
    .staticmethod("interaction_type_name")
    .def("interaction_type_short_name",&DotScorer::interaction_type_short_name)
    .staticmethod("interaction_type_short_name")
    .def("test", &DotScorer::test)
    ;

  // Export the vector indexing of objects that we'll use vectors for.
  // NOTE: Everything that is using scitbx::af::shared "flex" arrays is
  // automatically wrapped for us in ways that let them be used as standard
  // Python iterators so we don't need to add the wrapping.  We only need
  // to describe any std::vector values that we use.

  // Export the global functions
  def("closest_contact", closest_contact, "Point of closest contact and distance for dot on atom.");
  def("atom_charge", atom_charge, "Integer charge on an atom given its string charge description.");
  def("dot2srcCenter", dot2srcCenter, "Distance from dot to point on the source surface closest to target.");
  def("kissEdge2bullsEye", kissEdge2bullsEye, ".");
  def("annularDots", annularDots, ".");

  def("DotSpheres_test", DotSpheres_test, "Test all classes defined in DotSpheres.h.");
  def("SpatialQuery_test", SpatialQuery_test, "Test all classes defined in SpatialQuery.h.");
  def("Scoring_test", Scoring_test, "Test all classes defined in Scoring.h.");

  // Export the helper functions
  def("set_atom_i_seq", set_atom_i_seq, "Set the i_seq on an atom, required for Phantom Hydrogen processing.");
}
