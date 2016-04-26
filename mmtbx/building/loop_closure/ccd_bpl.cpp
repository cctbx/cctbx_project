#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>
#include <mmtbx/building/loop_closure/ccd.h>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_non_const_reference.hpp>
#include <boost/python/copy_const_reference.hpp>


namespace mmtbx { namespace building { namespace loop_closure {
namespace {

  void init_module()
  {
    using namespace boost::python;
    typedef return_value_policy<return_by_value> rbv;
    typedef return_value_policy<copy_const_reference> ccr;
    typedef return_value_policy<copy_non_const_reference> cncr;
    typedef default_call_policies dcp;

    class_<ccd_cpp> ("ccd_cpp", no_init)
      .def(init< scitbx::af::tiny<scitbx::vec3<double>, 3 >,
        iotbx::pdb::hierarchy::root & ,
        scitbx::af::tiny<size_t, 3>,
        const int&,
        const double&>(
          (arg("fixed_ref_atoms"),
           arg("moving_h"),
           arg("moving_ref_atoms_iseqs"),
           arg("max_number_of_iterations")=500,
           arg("needed_rmsd")=0.1)))
      .def("_find_angle", &ccd_cpp::_find_angle,
          (arg("axis_point_1"),
           arg("axis_point_2")))
      .def("_modify_angle", &ccd_cpp::_modify_angle,
          (arg("angle")))

      .add_property("early_exit",
          make_getter(&ccd_cpp::early_exit, rbv()),
          make_setter(&ccd_cpp::early_exit, dcp()))
      .add_property("resulting_rmsd",
          make_getter(&ccd_cpp::resulting_rmsd, rbv()),
          make_setter(&ccd_cpp::resulting_rmsd, dcp()))
      .add_property("needed_rmsd",
          make_getter(&ccd_cpp::needed_rmsd, rbv()),
          make_setter(&ccd_cpp::needed_rmsd, dcp()))
      .add_property("max_number_of_iterations",
          make_getter(&ccd_cpp::max_number_of_iterations, rbv()),
          make_setter(&ccd_cpp::max_number_of_iterations, dcp()))
      .add_property("r",
          make_getter(&ccd_cpp::r, rbv()),
          make_setter(&ccd_cpp::r, dcp()))
      .add_property("convergence_diff",
          make_getter(&ccd_cpp::convergence_diff, rbv()),
          make_setter(&ccd_cpp::convergence_diff, dcp()))
      .add_property("fixed_ref_atoms",
          make_getter(&ccd_cpp::fixed_ref_atoms, rbv()))
      .add_property("moving_ref_atoms_iseqs",
          make_getter(&ccd_cpp::moving_ref_atoms_iseqs, rbv()))
      .add_property("moving_h", make_getter(&ccd_cpp::moving_h, cncr()))
    ;
  }



} // namespace <anonymous>
}}} // namespace mmtbx::building::loop_closure

BOOST_PYTHON_MODULE(mmtbx_building_loop_closure_ext)
{
  mmtbx::building::loop_closure::init_module();
}
