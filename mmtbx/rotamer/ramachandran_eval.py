from __future__ import absolute_import, division, print_function

import boost_adaptbx.boost.python as bp
ext = bp.import_ext("mmtbx_validation_ramachandran_ext")
from mmtbx_validation_ramachandran_ext import rama_eval

# maps programatic name to file name
aminoAcids = {
    'general' : 'general',
    'glycine' : 'gly-sym',
    'proline' : 'pro',
    'prepro' : 'prepro',
}
aminoAcids_8000 = {
    'general' : 'general-noGPIVpreP',
    'glycine' : 'gly-sym',
    'cis-proline' : 'cispro',
    'trans-proline' : 'transpro',
    'pre-proline' : 'prepro-noGP',
    'isoleucine or valine' : 'ileval-nopreP',
}

#
# Why constants from ramalyze, i.e. res_types are not used here at all?
# They should be defined either here (preferably) and used everywhere or there.
#

class RamachandranEval:
  def __init__(self):
    self.rama_eval = rama_eval()

  def check_table_name(self, name):
    # This function take time to run. We have similar check in C++ and raising
    # Runtime error there if the name is not good.
    return name in aminoAcids_8000

  def evaluate(self, aaName, phiPsi):
    # assert self.check_table_name(aaName)
    return self.rama_eval.get_score(aaName, phiPsi[0],phiPsi[1])

  def evaluate_sites(self, aaName, phi_psi_i_seqs, sites_cart):
    # assert self.check_table_name(aaName)
    (phi, psi) = mmtbx.rotamer.phi_psi_from_sites(
      i_seqs=phi_psi_i_seqs,
      sites_cart=sites_cart)
    return self.evaluate(aaName, (phi,psi))
