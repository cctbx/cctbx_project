#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <iotbx/pdb/hybrid_36_c.h>
#include <string>
#include <stdexcept>

namespace {

  std::string
  hy36encode_wrapper(unsigned width, int value)
  {
    char result[16];
    const char* diag = hy36encode(width, value, result);
    if (diag) throw std::runtime_error(diag);
    return std::string(result);
  }

  int
  hy36decode_wrapper(unsigned width, std::string const& s)
  {
    int result;
    const char* diag = hy36decode(width, s.c_str(), s.size(), &result);
    if (diag) throw std::runtime_error(diag);
    return result;
  }

  unsigned
  hy36recode_width_4_all() // unit test
  {
    unsigned n_ok = 0;
    for(int value=-999;value<10000+2*26*36*36*36;value++) {
      char encoded[16];
      const char* diag = hy36encode(4U, value, encoded);
      if (!diag) {
        int decoded;
        const char* diag = hy36decode(4U, encoded, 4U, &decoded);
        if (!diag && decoded == value) {
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
