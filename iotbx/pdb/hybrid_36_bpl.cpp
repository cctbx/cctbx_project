#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <iotbx/pdb/hybrid_36_cpp.h>
#include <iotbx/pdb/hybrid_36_c.h>
#include <string>
#include <stdexcept>

namespace {

  std::string
  hy36encode_wrapper(unsigned width, int value)
  {
    return iotbx::pdb::hybrid_36::encode(width, value);
  }

  int
  hy36decode_wrapper(unsigned width, std::string const& s)
  {
    return iotbx::pdb::hybrid_36::decode(width, s);
  }

  unsigned
  hy36recode_width_4_all() // unit test
  {
    unsigned n_ok = 0;
    for(int value=-999;value<10000+2*26*36*36*36;value++) {
      char encoded[16];
      const char* errmsg = hy36encode(4U, value, encoded);
      if (!errmsg) {
        int decoded;
        const char* errmsg = hy36decode(4U, encoded, 4U, &decoded);
        if (!errmsg && decoded == value) {
          n_ok++;
        }
      }
    }
    return n_ok;
  }

} // namespace <anonymous>

namespace iotbx { namespace pdb { namespace boost_python {

  void
  wrap_hybrid_36()
  {
    using namespace boost::python;
    def("hy36encode", hy36encode_wrapper, (arg_("width"), arg_("value")));
    def("hy36decode", hy36decode_wrapper, (arg_("width"), arg_("s")));
    def("hy36recode_width_4_all", hy36recode_width_4_all);
  }

}}} // namespace iotbx::pdb::boost_python
