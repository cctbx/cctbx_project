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
  """
  Tables of side-chain chi angles of favored (or favored and allowed) rotamers.
  The files are large (up to 23 Mb, ~7 s to read all of them), so each table is
  read when first asked for, unless a list of residues is given, in which case
  those are read right away.
  """
  def __init__(self, rotamers, residues=None):
    assert rotamers in ["favored", "favored_allowed"]
    self.rotamers = rotamers
    self.rotamer_evaluator = RotamerEval()
    self.path = libtbx.env.find_in_repositories("chem_data/rotamer_chi_angles")
    self.result = {}
    self._flat = {}
    if(residues is not None):
      self.preload(resnames = residues)

  def _load_one(self, aa_code):
    if(aa_code in self.result): return
    try:
      chi_angles = easy_pickle.load(file_name = "/".join(
        [self.path, "%s_%s.pkl"%(aa_code, self.rotamers)]))
    except Exception:
      chi_angles = None
    self.result[aa_code] = [chi_angles]

  def preload(self, resnames):
    """
    Read tables for given residue names (those that are not amino acids are
    ignored). Use before passing this object to parallel workers, otherwise
    each of them reads the tables it needs.
    """
    requested = set([resname.strip().lower() for resname in resnames])
    for aa_code in aa_codes:
      if(aa_code in requested):
        self._load_one(aa_code = aa_code)
        self.get_chi_angles_flat(resname = aa_code)

  def get_chi_angles(self, resname):
    resname = resname.strip().lower()
    assert resname in aa_codes
    self._load_one(aa_code = resname)
    return self.result[resname][0]

  def get_chi_angles_flat(self, resname):
    """
    Same table as one array of numbers, sets of angles one after another,
    sorted by angles (first angle, then second, ...). Made once per residue
    type. This is what mmtbx_rotamer_fit_ext.fit takes without converting
    half a million Python lists for every residue; sets next to each other
    share leading angles, which the search uses not to repeat work.
    None if there is no table.
    """
    resname = resname.strip().lower()
    if(not resname in self._flat):
      chi_angles = self.get_chi_angles(resname = resname)
      if(chi_angles is None): self._flat[resname] = None
      else:
        from scitbx.array_family import flex
        if(len(chi_angles) == 0): self._flat[resname] = flex.double()
        else: # rows of equal length (flex checks)
          self._flat[resname] = flex.double(sorted(chi_angles)).as_1d()
    return self._flat[resname]

  def rotamer(self, residue):
    return self.rotamer_evaluator.evaluate_residue(residue)
