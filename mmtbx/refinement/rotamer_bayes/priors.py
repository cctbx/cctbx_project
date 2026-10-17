"""
Prior probabilities of side-chain rotamers, from the Top8000 distributions
in chem_data/rotarama_data (rota8000-*.data).

The tables hold, for each chi-angle cell, not the density but the fraction
of the distribution in cells of lower or equal density (the contour level
through the cell): 1 at the peak, 0.003 at the outlier contour, 0.02 at the
allowed contour. Sorting the cells by value, each step in value is the
probability mass of the cells at that level, so the masses are recovered by
differencing (ties share the step) and sum to the peak value, 1.

The cells are assigned to the rotamer wells of mmtbx/rotamer/
rotamer_names.props (first matching box, as mmtbx.rotamer.rotamer_eval.
RotamerID). Cells below the outlier contour form a single OUTLIER hypothesis
(mass 0.003 by definition) and cells above it in no box an UNCLASSIFIED one.
Each well's chi angles are the canonical ones of sidechain_angles.props if
listed, otherwise the table's peak cell.
"""
from __future__ import absolute_import, division, print_function
import os
import numpy as np
from mmtbx.rotamer import rotamer_eval
from mmtbx.rotamer.sidechain_angles import SidechainAngles

OUTLIER_THRESHOLD = 0.003
OUTLIER = "OUTLIER"
UNCLASSIFIED = "UNCLASSIFIED"

class well(object):
  def __init__(self, name, mass, chis):
    self.name = name
    self.mass = mass
    self.chis = chis   # degrees, as measureChiAngles (one per chi)

  def __repr__(self):
    return "well(%s, %.4f, %s)" % (self.name, self.mass,
      ",".join(["%.0f" % c for c in self.chis]))

def _read_table(file_name):
  """ [(cell centre tuple, value)] from a rota8000 .data file """
  a = np.loadtxt(file_name, comments="#", ndmin=2)
  return [(tuple(row[:-1]), row[-1]) for row in a]

def cell_masses(cells):
  """ Probability mass of each cell, by differencing the sorted contour
  levels (cells with equal values share the step). """
  by_value = {}
  for c, v in cells:
    by_value.setdefault(v, []).append(c)
  result = {}
  previous = 0.
  for v in sorted(by_value):
    cs = by_value[v]
    for c in cs:
      result[c] = (v - previous) / len(cs)
    previous = v
  return result

class rotamer_priors(object):
  """ Wells and their prior masses for each residue type with a table. """

  _cache = {}

  def __init__(self):
    self.rot_id = rotamer_eval.RotamerID()
    self.angles = SidechainAngles(False)
    self.data_dir = rotamer_eval.find_rotarama_data_dir()

  def wells(self, resname):
    """ {name: well} for resname (lower case; mse as met), including
    OUTLIER and UNCLASSIFIED; None if there is no table. """
    aa = resname.strip().lower()
    if(aa == "mse"): aa = "met"
    if(aa in self._cache): return self._cache[aa]
    table = rotamer_eval.aminoAcids.get(aa)
    if(table is None or aa not in self.rot_id.names):
      self._cache[aa] = None
      return None
    cells = _read_table(os.path.join(self.data_dir,
      "rota8000-%s.data" % table))
    masses = cell_masses(cells)
    coords = np.array([c for c, v in cells])
    values = np.array([v for c, v in cells])
    mass = np.array([masses[c] for c, v in cells])
    # first matching box, as RotamerID.identify (cell centres are wrapped)
    label = np.full(len(cells), -1)
    rots = self.rot_id.names[aa]
    for k, rot in enumerate(rots):
      b = np.array(rot.bounds, dtype=float).reshape(-1, 2)
      inside = np.all((coords >= b[:, 0]) & (coords <= b[:, 1]), axis=1)
      label[(label < 0) & inside] = k
    names = np.array([rots[k].rotamer_name if k >= 0 else UNCLASSIFIED
      for k in label], dtype=object)
    names[values < OUTLIER_THRESHOLD] = OUTLIER
    total = {}
    peak = {}
    for name in set(names.tolist()):
      sel = names == name
      total[name] = float(mass[sel].sum())
      i = np.nonzero(sel)[0][np.argmax(values[sel])]
      peak[name] = tuple(coords[i])
    result = {}
    for name, m in total.items():
      chis = None
      if(name not in (OUTLIER, UNCLASSIFIED)):
        chis = self.angles.get_rotamer_angles(aa, name)
      if(chis is None): chis = list(peak[name])
      result[name] = well(name=name, mass=m, chis=chis)
    self._cache[aa] = result
    return result

  def classify(self, resname, chis):
    """ Name of the well holding chis (degrees), OUTLIER or UNCLASSIFIED. """
    aa = resname.strip().lower()
    if(aa == "mse"): aa = "met"
    v = _evaluator().evaluate(aa, chis)
    if(v is None or v < OUTLIER_THRESHOLD): return OUTLIER
    wrapped = self.rot_id.wrap_chis(aa, chis, symmetry=True)
    return self.rot_id.identify(aa, wrapped) or UNCLASSIFIED

_ev = []
def _evaluator():
  if(len(_ev) == 0): _ev.append(rotamer_eval.RotamerEval())
  return _ev[0]
