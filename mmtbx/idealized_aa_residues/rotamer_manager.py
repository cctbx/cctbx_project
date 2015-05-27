from __future__ import division
from libtbx import easy_pickle
import libtbx.load_env
import math
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
  "mse", # XXX is ignored with rotamer named None
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
  def __init__(self):
    self.rotamer_evaluator = RotamerEval()
    path=libtbx.env.find_in_repositories("mmtbx/idealized_aa_residues/data")
    self.result = {}
    for aa_code in aa_codes:
      try:
        chi_angles = easy_pickle.load(
          file_name="/".join([path,"%s-coarse_step10.pickle"%aa_code]))
        # XXX Fix later
        chi_angles_ = []
        for chi in chi_angles:
          chi_ = []
          for ch in chi:
            chi_.append(ch*math.pi/180)
          chi_angles_.append(chi_)
        chi_angles = chi_angles_
        # XXX Fix later
      except Exception:
        chi_angles = None
      self.result.setdefault(aa_code,[]).extend([chi_angles])

  def get_chi_angles(self, resname):
    resname = resname.strip().lower()
    assert resname in aa_codes
    return self.result[resname][0]

  def rotamer(self, residue):
    return self.rotamer_evaluator.evaluate_residue(residue)
