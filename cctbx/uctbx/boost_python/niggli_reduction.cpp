#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <cctbx/uctbx/niggli_reduction.h>

namespace cctbx { namespace uctbx { namespace boost_python {
namespace {

  struct niggli_reduction_wrappers
  {
    typedef niggli_reduction w_t;

    // Return the 9 elements of r_inv as a Python tuple for use with
    // sgtbx.rot_mx(elems, den) in Python.
    static boost::python::tuple
    r_inv_as_tuple(w_t const& self)
    {
      scitbx::mat3<int> const& m = self.r_inv();
      return boost::python::make_tuple(
        m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8]);
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("niggli_reduction", no_init)
        .def(init<unit_cell const&,
                  optional<double, int> >((
          arg("unit_cell"),
          arg("relative_epsilon")=1.e-5,
          arg("iteration_limit")=1000)))
        .def("r_inv", &w_t::r_inv, ccr())
        .def("r_inv_as_tuple", r_inv_as_tuple)
        .def("n_iterations", &w_t::n_iterations)
        .def("as_unit_cell", &w_t::as_unit_cell)
      ;
    }
  };

} // namespace <anonymous>

  void wrap_niggli_reduction()
  {
    niggli_reduction_wrappers::wrap();
  }

}}} // namespace cctbx::uctbx::boost_python
