#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <mmtbx/alignment/alignment.h>
#include <mmtbx/alignment/align.h>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>


namespace mmtbx { namespace alignment {
namespace {

  struct pairwise_global_wrappers
  {
    typedef pairwise_global w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("pairwise_global", no_init)
        .def(init<std::string const&, std::string const&>((
          arg("seq1"), arg("seq2"))))
        .def_readonly("result1", &w_t::result1)
        .def_readonly("result2", &w_t::result2)
      ;
    }
  };

  struct align_global_wrappers
  {
    typedef align w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("align", no_init)
        .def(init<std::string const&,
                  std::string const&,
                  std::string const&,
              float, float, std::string const&>((
          arg("seq_a"),
          arg("seq_b"),
          arg("style")="global",
          arg("gap_opening_penalty")=1,
          arg("gap_extension_penalty")=1,
          arg("similarity_function")="identity")))
        .def(init<af::const_ref<std::string> const&,
                  af::const_ref<std::string> const&,
                  std::string const&,
                  float, float, std::string const&>((
          arg("seq_a"),
          arg("seq_b"),
          arg("style")="global",
          arg("gap_opening_penalty")=1,
          arg("gap_extension_penalty")=1,
          arg("similarity_function")="identity")))
        .add_property("M", make_getter(&w_t::M, rbv()))
        .add_property("D", make_getter(&w_t::D, rbv()))
        .add_property("I", make_getter(&w_t::I, rbv()))
        .add_property("E", make_getter(&w_t::E, rbv()))
      ;
    }
  };


  void init_module()
  {
    pairwise_global_wrappers::wrap();
    align_global_wrappers::wrap();
  }

} // namespace <anonymous>
}} // namespace mmtbx::alignment

BOOST_PYTHON_MODULE(mmtbx_alignment_ext)
{
  mmtbx::alignment::init_module();
}
