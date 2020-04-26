from __future__ import absolute_import, division, print_function
import boost_adaptbx.python
ext = boost_adaptbx.python.import_ext("cctbx_eltbx_neutron_ext")
from cctbx_eltbx_neutron_ext import *

boost_adaptbx.python.inject(ext.neutron_news_1992_table_iterator, boost_adaptbx.python.py3_make_iterator)
