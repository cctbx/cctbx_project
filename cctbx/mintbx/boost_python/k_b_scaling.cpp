/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Nov: Created (Ralf W. Grosse-Kunstleve)
 */

#include <cctbx/mintbx/k_b_scaling.h>
#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>

namespace cctbx { namespace mintbx { namespace boost_python {
namespace {

  struct k_b_scaling_target_and_gradients_wrappers
  {
    typedef k_b_scaling_target_and_gradients<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("k_b_scaling_target_and_gradients", no_init)
        .def(init<uctbx::unit_cell const&,
                  af::const_ref<miller::index<> > const&,
                  af::const_ref<int> const&,
                  af::const_ref<double> const&,
                  af::const_ref<double> const&,
                  double,
                  double,
                  bool,
                  bool>())
        .def(init<uctbx::unit_cell const&,
                  af::const_ref<miller::index<> > const&,
                  af::const_ref<int> const&,
                  af::const_ref<double> const&,
                  af::const_ref<double> const&,
                  double,
                  scitbx::sym_mat3<double> const&,
                  bool,
                  bool>())
        .def("target", &w_t::target)
        .def("gradient_k", &w_t::gradient_k)
        .def("anisotropic_flag", &w_t::anisotropic_flag)
        .def("gradient_b_iso", &w_t::gradient_b_iso)
        .def("gradients_b_cif", &w_t::gradients_b_cif, ccr())
      ;
    }
  };

} // namespace <anonymous>

  void wrap_k_b_scaling()
  {
    k_b_scaling_target_and_gradients_wrappers::wrap();
  }

}}} // namespace cctbx::mintbx::boost_python
