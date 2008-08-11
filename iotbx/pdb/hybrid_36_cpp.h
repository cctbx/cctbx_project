/* If you change the include guards, please be sure to also change the
   namespaces below. Otherwise your project will clash with the original
   iotbx declarations and definitions.
 */
#ifndef IOTBX_PDB_HYBRID_36_CPP_H
#define IOTBX_PDB_HYBRID_36_CPP_H

#include <string>

namespace iotbx { namespace pdb {

//! Hybrid-36 C++ wrappers.
/*!
    See also: http://cci.lbl.gov/hybrid_36/
 */
namespace hybrid_36 {

  std::string
  encode(unsigned width, int value);

  int
  decode(unsigned width, const char* s, unsigned s_size);

  int
  decode(unsigned width, std::string const& s);

}}} // namespace iotbx::pdb::hybrid_36

#endif // IOTBX_PDB_HYBRID_36_CPP_H
