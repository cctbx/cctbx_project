#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <iotbx/pdb/hierarchy.h>

namespace iotbx { namespace pdb {
namespace {

  struct atom_wrappers
  {
    typedef atom w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("atom", no_init)
      ;
    }
  };

  void
  wrap_all()
  {
    using namespace boost::python;
    atom_wrappers::wrap();
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_hierarchy() { wrap_all(); }

}}} // namespace iotbx::pdb::boost_python
