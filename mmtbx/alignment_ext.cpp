#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <mmtbx/alignment.h>

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

  void init_module()
  {
    pairwise_global_wrappers::wrap();
  }

} // namespace <anonymous>
}} // namespace mmtbx::alignment

BOOST_PYTHON_MODULE(mmtbx_alignment_ext)
{
  mmtbx::alignment::init_module();
}
