from scitbx.source_generators.array_family import generate_all
from libtbx.utils import warn_if_unexpected_md5_hexdigest
from scitbx.source_generators import flex_fwd_h
import os

if (self.env.is_ready_for_build()):
  message_template = '  Generating C++ files in:\n    "%s"'

  # array_family
  target_dir = self.env.under_build("include/scitbx/array_family/detail")
  if (not os.path.isdir(target_dir)):
    os.makedirs(target_dir)
  generate_all.refresh(array_family=os.path.dirname(target_dir))
  #
  warn_if_unexpected_md5_hexdigest(
    path=self.env.under_dist(
      module_name="boost", path="boost/random/mersenne_twister.hpp"),
    expected_md5_hexdigests=[
      "a2533c79a21f0f773f2e0d29a37371b0", # CVS revision 1.21
      "378432b5b280c9c0c894f7c80f0dad92", # CVS revision 1.20
      "1fe430a94f330e36dca14ef2a553dba1", # SVN revision 37030 (== CVS 1.20)
    ],
    hints=[
      "  Files to review:",
      "    scitbx/include/scitbx/random.h",
      "    scitbx/libtbx_refresh.py"])

  # flex_fwd.h
  target_dir = self.env.under_build("include/scitbx/array_family/boost_python")
  print message_template % target_dir
  if not os.path.isdir(target_dir):
    os.makedirs(target_dir)
  flex_fwd_h.run(target_dir)
