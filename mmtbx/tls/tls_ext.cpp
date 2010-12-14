#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <mmtbx/tls/tls.h>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <scitbx/boost_python/is_polymorphic_workaround.h>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python.hpp>

SCITBX_BOOST_IS_POLYMORPHIC_WORKAROUND(mmtbx::tls::common)

namespace mmtbx { namespace tls {
  namespace bp = boost::python;

namespace {
  boost::python::tuple
  getinitargs(tlso<> const& self)
  {
    return boost::python::make_tuple(self.t, self.l, self.s);
  }

  void init_module()
  {
    using namespace boost::python;
    using boost::python::arg;
    typedef return_value_policy<return_by_value> rbv;

    class_<uaniso_from_tls>("uaniso_from_tls",
                             init<sym_mat3<double> const&,
                                  sym_mat3<double> const&,
                                  mat3<double> const&,
                                  vec3<double> const&,
                                  vec3<double> const&>())
      .def("u", &uaniso_from_tls::u)
    ;

    class_<tls_ls_derivative_coefficients>("tls_ls_derivative_coefficients",
                           init<vec3<double> const&,
                                af::shared<vec3<double> > const&,
                                af::shared<double> const&>(
                                   (arg("origin"),
                                    arg("sites_cart"),
                                    arg("u_iso"))))
      .add_property("a", make_getter(&tls_ls_derivative_coefficients::a, rbv()))
      .add_property("b", make_getter(&tls_ls_derivative_coefficients::b, rbv()))
    ;

    class_<d_target_d_tls>("d_target_d_tls",
                           init<af::shared<vec3<double> > const&,
                                vec3<double> const&,
                                af::shared<sym_mat3<double> > const&,
                                bool,
                                bool>(
                                   (arg("sites"),
                                    arg("origin"),
                                    arg("d_target_d_uaniso"),
                                    arg("scale_l_and_s"),
                                    arg("use_trace_s_zero_constraint"))))
      .def("grad_T", &d_target_d_tls::grad_T)
      .def("grad_L", &d_target_d_tls::grad_L)
      .def("grad_S", &d_target_d_tls::grad_S)
    ;

    class_<tls_from_uaniso_target_and_grads>("tls_from_uaniso_target_and_grads",
                             init<sym_mat3<double> const&,
                                  sym_mat3<double> const&,
                                  mat3<double> const&,
                                  vec3<double> const&,
                                  af::shared<vec3<double> > const&,
                                  af::shared<sym_mat3<double> > const&>())
      .def("target", &tls_from_uaniso_target_and_grads::target)
      .def("grad_T", &tls_from_uaniso_target_and_grads::grad_T)
      .def("grad_L", &tls_from_uaniso_target_and_grads::grad_L)
      .def("grad_S", &tls_from_uaniso_target_and_grads::grad_S)
    ;
    class_<tls_from_uiso_target_and_grads>("tls_from_uiso_target_and_grads",
                             init<double const&,
                                  sym_mat3<double> const&,
                                  vec3<double> const&,
                                  vec3<double> const&,
                                  af::shared<vec3<double> > const&,
                                  af::shared<double> const&>())
      .def("target", &tls_from_uiso_target_and_grads::target)
      .def("grad_T", &tls_from_uiso_target_and_grads::grad_T)
      .def("grad_L", &tls_from_uiso_target_and_grads::grad_L)
      .def("grad_S", &tls_from_uiso_target_and_grads::grad_S)
    ;

    typedef return_value_policy<return_by_value> rbv;
    class_<tlso<> >("tlso")
      .def(init<scitbx::sym_mat3<double> const&,
                scitbx::sym_mat3<double> const&,
                scitbx::mat3<double> const&,
                scitbx::vec3<double> const& >((arg("t"),arg("l"),arg("s"),
                                               arg("origin"))))
      .add_property("t",      make_getter(&tlso<>::t,      rbv()))
      .add_property("l",      make_getter(&tlso<>::l,      rbv()))
      .add_property("s",      make_getter(&tlso<>::s,      rbv()))
      .add_property("origin", make_getter(&tlso<>::origin, rbv()))
      .enable_pickling()
      .def("__getinitargs__", getinitargs)
    ;
    def("uaniso_from_tls_one_group",
         (af::shared<sym_mat3<double> >(*)
               (tlso<double>,
                af::shared<vec3<double> > const&)) uaniso_from_tls_one_group,
                                                          (arg("tlso"),
                                                           arg("sites_cart")))
   ;

   class_<tls_parts_one_group>("tls_parts_one_group",
                             init<tlso<double>,
                                  af::shared<vec3<double> > const&>())
      .def("ala",    &tls_parts_one_group::ala)
      .def("assa",   &tls_parts_one_group::assa)
      .def("u_cart", &tls_parts_one_group::u_cart)
      .def("t",      &tls_parts_one_group::t)
      .def("r",      &tls_parts_one_group::r)
    ;
   class_<tls_parts_one_group_as_b_iso>("tls_parts_one_group_as_b_iso",
                             init<tlso<double>,
                                  af::shared<vec3<double> > const&>())
      .def("ala",    &tls_parts_one_group_as_b_iso::ala)
      .def("assa",   &tls_parts_one_group_as_b_iso::assa)
      .def("b_iso", &tls_parts_one_group_as_b_iso::b_iso)
      .def("t",      &tls_parts_one_group_as_b_iso::t)
    ;

   class_<common>("common",init<sym_mat3<double> const&,
                                sym_mat3<double> const&,
                                optional<double> >())
      .def("t", &common::t)
      .def("branch_0",       &common::get_branch_0)
      .def("branch_1",       &common::get_branch_1)
      .def("branch_1_1",     &common::get_branch_1_1)
      .def("branch_1_2",     &common::get_branch_1_2)
      .def("branch_1_2_1",   &common::get_branch_1_2_1)
      .def("branch_1_2_2",   &common::get_branch_1_2_2)
      .def("branch_1_2_3",   &common::get_branch_1_2_3)
      .def("branch_1_2_3_1", &common::get_branch_1_2_3_1)
      .def("branch_1_2_3_2", &common::get_branch_1_2_3_2)
   ;

   def("ls_target_from_iso_tls", (double(*)(
     double const&,
     sym_mat3<double> const&,
     vec3<double> const&,
     vec3<double> const&,
     af::shared<vec3<double> > const&,
     af::shared<double> const&)) ls_target_from_iso_tls, (arg("t"),arg("l"),
       arg("s"),arg("origin"),arg("sites_cart"),arg("u_isos")))
   ;

   def("t_from_u_cart", (sym_mat3<double>(*)(af::shared<sym_mat3<double> > const&, double)) t_from_u_cart, (arg("u_cart"),arg("small")))
   ;
   def("t_from_u_cart", (sym_mat3<double>(*)(af::shared<double> const&, double))            t_from_u_cart, (arg("u_iso"),arg("small")))
   ;

  }

} // namespace <anonymous>
}} // namespace mmtbx::tls

BOOST_PYTHON_MODULE(mmtbx_tls_ext)
{
  mmtbx::tls::init_module();
}
