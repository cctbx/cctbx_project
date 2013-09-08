#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include <mmtbx/secondary_structure/dssp.hpp>

namespace mmtbx { namespace secondary_structure { namespace dssp {
namespace {
  void wrap_dssp_hbond ()
  {
    using namespace boost::python;
    def("get_o_n_hbond_energy", get_o_n_hbond_energy, (
      arg("O"), arg("N"), arg("nh_bond_length")=1.01));
  }
}

namespace boost_python {
  void wrap_dssp ()
  {
    wrap_dssp_hbond();
  }
}
}}}

BOOST_PYTHON_MODULE(mmtbx_dssp_ext)
{
  mmtbx::secondary_structure::dssp::boost_python::wrap_dssp();
}
