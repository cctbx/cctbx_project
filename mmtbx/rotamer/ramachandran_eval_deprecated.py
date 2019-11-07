from __future__ import absolute_import, division, print_function
import mmtbx.rotamer
from mmtbx.rotamer.n_dim_table import NDimTable
from mmtbx.rotamer.rotamer_eval import find_rotarama_data_dir
from mmtbx.rotamer.rotamer_eval import open_rotarama_dlite
from libtbx import easy_pickle
from libtbx.utils import Sorry
import weakref
import os

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

class RamachandranEval:

  # This is shared among all instances of RamachandranEval -- a class variable.
  # It holds a LOT of read-only data, so this helps save memory.
  aaTables = {} # maps "his" to a NDimTable object for histidine, etc.

  def __init__(self):
    main_aaTables = RamachandranEval.aaTables
    self.aaTables = {}
    for aa,ndt_weakref in main_aaTables.items():
      # convert existing weak references to strong references
      self.aaTables[aa] = ndt_weakref()
    rama_data_dir = find_rotarama_data_dir()
    target_db = open_rotarama_dlite(rotarama_data_dir=rama_data_dir)
    no_update = os.path.exists(os.path.join(rama_data_dir, "NO_UPDATE"))
    for aa, aafile in aminoAcids_8000.items():
      if (self.aaTables.get(aa) is not None): continue
      data_file = "rama8000-"+aafile+".data"
      pickle_file = "rama8000-"+aafile+".pickle"
      pair_info = target_db.pair_info(
        source_path=data_file,
        target_path=pickle_file,
        path_prefix=rama_data_dir)
      if (((pair_info.needs_update) and (not no_update)) or not
          os.path.exists(os.path.join(rama_data_dir, pickle_file))):
        raise Sorry(
          "chem_data/rotarama_data/*.pickle files are missing or out of date.\n"
          "  Please run\n"
          "    mmtbx.rebuild_rotarama_cache\n"
          "  to resolve this problem.\n")
      ndt = easy_pickle.load(file_name=os.path.join(
        rama_data_dir, pair_info.target.path))
      self.aaTables[aa] = ndt
      main_aaTables[aa] = weakref.ref(ndt)

  def evaluate(self, aaName, phiPsi):
    '''Evaluates the protein backbone conformation from 0.0 (worst) to 1.0 (best).

    Values below 0.0005 are considered outliers for the general case,
    or values below 0.002 for glycine, proline, and pre-proline.
    The field aaName should be one of "general", "glycine", "proline", or "prepro".
    If the aaName is not recognized, returns None.
    phiPsi is a list or tuple of angles in degrees(?): [phi, psi]'''
    ndt = self.aaTables.get(aaName.lower())
    if (ndt is None): return None
    return ndt.valueAt(phiPsi)

  def evaluate_sites(self, aaName, phi_psi_i_seqs, sites_cart):
    assert (aaName in ["general",
                       "glycine",
                       "cis-proline",
                       "trans-proline",
                       "pre-proline",
                       "isoleucine or valine"])
    (phi, psi) = mmtbx.rotamer.phi_psi_from_sites(
      i_seqs=phi_psi_i_seqs,
      sites_cart=sites_cart)
    return self.evaluate(aaName, (phi,psi))
