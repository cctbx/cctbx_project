from __future__ import absolute_import, division, print_function
from libtbx import easy_pickle
import libtbx.load_env
from mmtbx.rotamer.rotamer_eval import RotamerEval

aa_codes = [
  "ala",
  "asn",
  "asp",
  "cys",
  "gln",
  "glu",
  "gly",
  "his",
  "ile",
  "leu",
  "met",
  "mse",
  "phe",
  "pro", # XXX BAD all-rotamers files
  "ser",
  "thr",
  "trp",
  "tyr",
  "val",
  "arg",
  "lys"
]

class load(object):
  def __init__(self, rotamers, residues=None):
    if(residues is not None):
      for i, residue in enumerate(residues):
        residues[i] = residue.lower()
    assert rotamers in ["favored", "favored_allowed"]
    self.rotamer_evaluator = RotamerEval()
    path=libtbx.env.find_in_repositories("chem_data/rotamer_chi_angles")
    self.result = {}
    for aa_code in aa_codes:
      if(residues is not None):
        if(not aa_code in residues): continue
      try:
        if(rotamers=="favored"):
          f = "%s_favored.pkl"%aa_code
        elif(rotamers=="favored_allowed"):
          f = "%s_favored_allowed.pkl"%aa_code
        else: raise RuntimeError("Not implemented: rotamers=%s"%rotamers)
        chi_angles = easy_pickle.load(file_name = "/".join([path,f]))
      except Exception:
        chi_angles = None
      self.result.setdefault(aa_code,[]).extend([chi_angles])

  def get_chi_angles(self, resname):
    resname = resname.strip().lower()
    assert resname in aa_codes
    return self.result[resname][0]

  def rotamer(self, residue):
    return self.rotamer_evaluator.evaluate_residue(residue)
