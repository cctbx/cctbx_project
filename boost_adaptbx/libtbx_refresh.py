from __future__ import division
import re

# Find Boost version at configuration time by parsing boost/version.hpp
try:
  boost_root = self.env.find_in_repositories('boost',
                                             return_relocatable_path=True)
  version_pat = re.compile(r'^ \s* \#define \s+ BOOST_VERSION \s+ (\d+) \s* $',
                           re.X)
  self.env.boost_version = -1
  with (boost_root / 'boost' / 'version.hpp').open() as lines:
    for li in lines:
      m = version_pat.search(li)
      if m:
        self.env.boost_version = int(m.group(1))
except IOError:
  pass
