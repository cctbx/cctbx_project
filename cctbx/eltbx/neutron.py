from __future__ import absolute_import, division, print_function
import boost.python
ext = boost.python.import_ext("cctbx_eltbx_neutron_ext")
from cctbx_eltbx_neutron_ext import *

boost.python.inject(ext.neutron_news_1992_table_iterator, boost.python.py3_make_iterator)
