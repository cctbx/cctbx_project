from __future__ import absolute_import, division, print_function
from scitbx.source_generators.flex_fwd_h import write

this = "cctbx.source_generators.flex_fwd_h"

common_code = """\
#ifndef SMTBX_BOOST_PYTHON_FLEX_FWD_H
#define SMTBX_BOOST_PYTHON_FLEX_FWD_H

#include <cctbx/boost_python/flex_fwd.h>

%s#endif // SMTBX_BOOST_PYTHON_FLEX_FWD_H
"""

full_code = """\
%s

#if defined(__sgi) && !defined(__GNUC__)

namespace smtbx { namespace boost_python {

  template <typename T>
  struct flex_fwd
  {
  };

  inline void
  flex_fwd_types()
  {
  }

}} // namespace cctbx::boost_python

#endif // defined(__sgi) && !defined(__GNUC__)

"""

def run(target_dir):
  write(this, target_dir, common_code, full_code)

if (__name__ == "__main__"):
  run(".")
