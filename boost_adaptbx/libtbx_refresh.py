from __future__ import absolute_import, division, print_function
import os
import re
import sys

# Find Boost version at configuration time by parsing boost/version.hpp
try:
  self.env.boost_version = -1
  boost_root = self.env.find_in_repositories('boost')
  # boost directory in cctbx_project is picked up by find_in_repositories
  if not os.path.isdir(boost_root) or 'cctbx_project' in boost_root:
    boost_root = os.path.join(sys.prefix, 'include')
    if sys.platform == 'darwin' and 'python.app' in boost_root:
      boost_root = os.path.join(boost_root.split('python.app')[0], 'include')
    elif sys.platform == 'win32':
      boost_root = os.path.join(sys.prefix, 'Library', 'include')
  version_pat = re.compile(r'^ \s* \#define \s+ BOOST_VERSION \s+ (\d+) \s* $',
                            re.X)
  with open(os.path.join(boost_root, 'boost', 'version.hpp')) as lines:
    for li in lines:
      m = version_pat.search(li)
      if m:
        self.env.boost_version = int(m.group(1))
except IOError:
  pass
except TypeError:
  if boost_root is None: pass
