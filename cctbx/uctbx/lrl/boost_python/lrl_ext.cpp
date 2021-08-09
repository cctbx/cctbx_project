#include <boost/python/operators.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/overloads.hpp>

#include "LRL_LatticeMatcher.h"
#include "LRL_ReadLatticeData.h"
#include "LRL_Cell.h"
#include "MatS6.h"
#include "LatticeConverter.h"
#include "S6.h"
#include "DC.h"

namespace cctbx { namespace uctbx { namespace lrl { namespace boost_python {

using namespace boost::python;
typedef return_internal_reference<> rir;


static void wrap()
{
  class_<LRL_ReadLatticeData>("LRL_ReadLatticeData")
    .def <void (LRL_ReadLatticeData::*)(const std::string&)>
      ("CellReader", &LRL_ReadLatticeData::CellReader)
    .def("GetLattice", &LRL_ReadLatticeData::GetLattice)
    .def("GetCell", &LRL_ReadLatticeData::GetCell)
  ;

  class_<LRL_Cell>("LRL_Cell")
    .def(self_ns::str(self_ns::self))
  ;

  class_<LRL_Cell_Degrees>("LRL_Cell_Degrees")
    .def(init<LRL_Cell>())
    .def(init<S6>()) // this works via an implicit cast to LRL_Cell
    .def(self_ns::str(self_ns::self))
  ;

  class_<LRL_LatticeMatcher>("LRL_LatticeMatcher")
    // Here we wrap overloaded instance methods (contrast to SellingReduceCell
    // below); thus we static_cast to a member pointer
    .def(
        "SetReferenceLattice",
        static_cast<void (LRL_LatticeMatcher::*)(const S6&)>(
            &LRL_LatticeMatcher::SetReferenceLattice
        )
    )
    .def(
        "MatchReference",
        static_cast<S6 (LRL_LatticeMatcher::*)(const S6&) const>(
            &LRL_LatticeMatcher::MatchReference
        )
    )
    .def(
        "MatchReference",
        static_cast<
          std::vector<S6> (LRL_LatticeMatcher::*)(const std::vector<S6>&) const
        >(&LRL_LatticeMatcher::MatchReference)
    )
    .def("GetMaxRadius", &LRL_LatticeMatcher::GetMaxRadius)
  ;

  class_<MatS6>("MatS6")
    .def(self_ns::str(self_ns::self))
    .def("Inverse", &MatS6::Inverse)
    .def(self_ns::self * other<S6>())
    .staticmethod("Inverse")
  ;

  class_<S6>("S6")
    .def(init<LRL_Cell>())
    .def(self_ns::str(self_ns::self))
  ;


  /*
  This wraps an overloaded static member function by static_casting
  the function pointer to the correct call signature
  */
  class_<LatticeConverter>("LatticeConverter")
    .def(
        "SellingReduceCell",
        static_cast<LRL_Cell (*)(const std::string&, const LRL_Cell&, MatS6&)>(
            &LatticeConverter::SellingReduceCell
        )
    )
    .def(
        "SellingReduceCell",
        static_cast<LRL_Cell (*)(const std::string&, const LRL_Cell&)>(
            &LatticeConverter::SellingReduceCell
        )
    )
    .staticmethod("SellingReduceCell")
  ;

  class_<DC>("DC")
    .def(init<S6>())
    .def("DistanceBetween", &DC::DistanceBetween)
    .staticmethod("DistanceBetween")
  ;


}



}}}} // namespace cctbx::uctbx::lrl::boost_python


BOOST_PYTHON_MODULE(cctbx_uctbx_lrl_ext)
{
  using namespace cctbx::uctbx::lrl::boost_python;
  using namespace boost::python;

  wrap();



}
