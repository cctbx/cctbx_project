#include <boost/python/class.hpp>
#include <cctbx/sgtbx/find_affine.h>

namespace cctbx { namespace sgtbx {

namespace {

  struct find_affine_wrappers
  {
    typedef find_affine w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("find_affine", no_init)
        .def(init<space_group const&, optional<int> >())
        .def("cb_mx", &w_t::cb_mx)
      ;
    }
  };

} // namespace <anoymous>

namespace boost_python {

  void wrap_find_affine()
  {
    find_affine_wrappers::wrap();
  }

}}} // namespace cctbx::sgtbx::boost_python
