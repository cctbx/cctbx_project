from __future__ import absolute_import, division, print_function
import re

# Find Boost version at configuration time by parsing boost/version.hpp
try:
  self.env.boost_version = -1
  boost_root = self.env.find_in_repositories('boost',
                                             return_relocatable_path=True)
  version_pat = re.compile(r'^ \s* \#define \s+ BOOST_VERSION \s+ (\d+) \s* $',
                           re.X)
  with (boost_root / 'boost' / 'version.hpp').open() as lines:
    for li in lines:
      m = version_pat.search(li)
      if m:
        self.env.boost_version = int(m.group(1))
except IOError:
  pass
except TypeError:
  if boost_root is None: pass
