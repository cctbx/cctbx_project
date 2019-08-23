
#include <simtbx/nanoBragg/nanoBragg.h>

namespace simtbx {
namespace nanoBragg {

class diffBragg: public nanoBragg{
  public:
  inline
  diffBragg(const dxtbx::model::Detector& detector, const dxtbx::model::Beam& beam,
            int verbose, int panel_id = 0):
            nanoBragg(detector, beam, verbose, panel_id)
            {}
};
}}
