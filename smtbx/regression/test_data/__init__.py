from __future__ import division
import os

test_data_dir = os.path.dirname(__file__)

class fnames(object):
  thpp_ins = os.path.join(test_data_dir, 'thpp.ins')
  thpp_hkl = os.path.join(test_data_dir, 'thpp.hkl')
  thpp_cif = os.path.join(test_data_dir, 'thpp.cif')
  thpp_out = os.path.join(test_data_dir, 'thpp_out.cif')

  sucrose_p1_res = os.path.join(test_data_dir, 'sucrose_p1.res')
