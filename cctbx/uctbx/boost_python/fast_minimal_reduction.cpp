#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <cctbx/uctbx/fast_minimal_reduction.h>

namespace cctbx { namespace uctbx { namespace boost_python {

namespace {

  struct fast_minimal_reduction_wrappers
  {
    typedef fast_minimal_reduction<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("fast_minimal_reduction", no_init)
        .def(init<unit_cell const&, optional<std::size_t, std::size_t> >())
        .def("as_gruber_matrix", &w_t::as_gruber_matrix)
        .def("as_niggli_matrix", &w_t::as_niggli_matrix)
        .def("as_sym_mat3", &w_t::as_sym_mat3)
        .def("as_unit_cell", &w_t::as_unit_cell)
        .def("expected_cycle_limit", &w_t::expected_cycle_limit)
        .def("iteration_limit", &w_t::iteration_limit)
        .def("r_inv", &w_t::r_inv, ccr())
        .def("n_iterations", &w_t::n_iterations)
        .def("had_expected_cycle", &w_t::had_expected_cycle)
        .def("type", &w_t::type)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_fast_minimal_reduction()
  {
    fast_minimal_reduction_wrappers::wrap();
  }

}}} // namespace cctbx::sgtbx::boost_python
