#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <cctbx/uctbx/fast_minimum_reduction.h>

namespace cctbx { namespace uctbx { namespace boost_python {

namespace {

  struct fast_minimum_reduction_wrappers
  {
    typedef fast_minimum_reduction<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("fast_minimum_reduction", no_init)
        .def(init<unit_cell const&,
                  optional<std::size_t, double, std::size_t> >())
        .def("as_gruber_matrix", &w_t::as_gruber_matrix)
        .def("as_niggli_matrix", &w_t::as_niggli_matrix)
        .def("as_sym_mat3", &w_t::as_sym_mat3)
        .def("as_unit_cell", &w_t::as_unit_cell)
        .def("iteration_limit", &w_t::iteration_limit)
        .def("multiplier_significant_change_test",
             &w_t::multiplier_significant_change_test)
        .def("min_n_no_significant_change", &w_t::min_n_no_significant_change)
        .def("r_inv", &w_t::r_inv, ccr())
        .def("n_iterations", &w_t::n_iterations)
        .def("termination_due_to_significant_change_test",
             &w_t::termination_due_to_significant_change_test)
        .def("type", &w_t::type)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_fast_minimum_reduction()
  {
    fast_minimum_reduction_wrappers::wrap();
  }

}}} // namespace cctbx::sgtbx::boost_python
