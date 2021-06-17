from __future__ import absolute_import, division, print_function
import boost_adaptbx.boost.python as bp
ext = bp.import_ext("cctbx_eltbx_neutron_ext")
from cctbx_eltbx_neutron_ext import *

bp.inject(ext.neutron_news_1992_table_iterator, bp.py3_make_iterator)
