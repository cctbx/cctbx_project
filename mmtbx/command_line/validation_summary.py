
"""
Convenience tool for collecting validation statistics with minimal overhead.
"""

from __future__ import absolute_import, division, print_function
from mmtbx.validation import molprobity
import iotbx.pdb
from libtbx import slots_getstate_setstate, Auto
from libtbx.utils import Sorry, Usage
from libtbx import str_utils
from libtbx import easy_mp
import os
import sys
from six.moves import range
import mmtbx.model

def summary(pdb_file=None,
            pdb_hierarchy=None,
            crystal_symmetry=None):
  header_info = None
  if (pdb_hierarchy is None):
    assert (pdb_file is not None)
    pdb_in = iotbx.pdb.input(pdb_file)
    pdb_hierarchy = pdb_in.construct_hierarchy()
    pdb_hierarchy.atoms().reset_i_seq()
    header_info = molprobity.pdb_header_info(
      pdb_file=pdb_file)
    crystal_symmetry=pdb_in.crystal_symmetry()
  else :
    assert (pdb_file is None)
  #
  assert crystal_symmetry is not None

  cache = pdb_hierarchy.atom_selection_cache()
  sel = cache.selection('protein')
  pdb_hierarchy = pdb_hierarchy.select(sel)
  #
  model = pdb_hierarchy.as_model_manager(crystal_symmetry = crystal_symmetry)
  return molprobity.molprobity(
    model=model,
    keep_hydrogens=False,
    header_info=header_info).summarize()

class parallel_driver(object):
  """
  Simple wrapper for passing to easy_mp.pool_map.
  """
  def __init__(self, pdb_hierarchy, crystal_symmetry):
    self.pdb_hierarchy = pdb_hierarchy
    self.crystal_symmetry = crystal_symmetry

  def __call__(self, i_model):
    import iotbx.pdb.hierarchy
    model_hierarchy = iotbx.pdb.hierarchy.root()
    model = self.pdb_hierarchy.models()[i_model].detached_copy()
    model.id = ""
    model_hierarchy.append_model(model)
    return summary(pdb_hierarchy=model_hierarchy,
        crystal_symmetry=self.crystal_symmetry)

molprobity_stat_labels = [
  "Ramachandran Outliers",
  "Ramachandran Favored",
  "Rotamer Outliers",
  "C-beta Outliers",
  "Clashscore",
  "MolProbity Score",
]

class ensemble(slots_getstate_setstate):
  """
  MolProbity validation results for an ensemble of models.  Note that the
  number of atoms in each model is not necessarily consistent.
  """

  __slots__ = [
    "rama_outliers",
    "rama_favored",
    "rotamer_outliers",
    "c_beta_deviations",
    "clashscore",
    "mpscore",
  ]

  def __init__(self, pdb_hierarchy, n_models, crystal_symmetry, nproc=Auto):
    assert (len(pdb_hierarchy.models()) == n_models)
    validate = parallel_driver(pdb_hierarchy, crystal_symmetry)
    summaries = easy_mp.pool_map(
      processes=nproc,
      fixed_func=validate,
      args=range(n_models))
    for name in self.__slots__ :
      array = []
      for s in summaries :
        array.append(getattr(s, name))
      setattr(self, name, array)

  def show(self, out=None, prefix="", show_percentiles=None):
    if (out is None):
      out = sys.stdout
    def min_max_mean(array):
      if (len(array) == 0) or (array.count(None) == len(array)):
        return (None, None, None)
      else :
        return min(array), max(array), sum(array) / len(array)
    def fs(format, value):
      return str_utils.format_value(format, value, replace_none_with=("(none)"))
    def format_all(format, array):
      min, max, mean = min_max_mean(array)
      return "%s %s %s" % (fs(format, min), fs(format, max), fs(format, mean))
    print("%s                           min    max   mean" % prefix, file=out)
    print("%sRamachandran outliers = %s %%" % (prefix,
      format_all("%6.2f", self.rama_outliers)), file=out)
    print("%s             favored  = %s %%" % (prefix,
      format_all("%6.2f", self.rama_favored)), file=out)
    print("%sRotamer outliers      = %s %%" % (prefix,
      format_all("%6.2f", self.rotamer_outliers)), file=out)
    print("%sC-beta deviations     = %s" % (prefix,
      format_all("%6d", self.c_beta_deviations)), file=out)
    print("%sClashscore            = %s" % (prefix,
      format_all("%6.2f", self.clashscore)), file=out)
    if (self.mpscore is not None):
      print("%sMolprobity score      = %s" % (prefix,
        format_all("%6.2f", self.mpscore)), file=out)

def run(args, out=sys.stdout):
  import optparse
  if (len(args) == 0) or ("--help" in args):
    raise Usage("""
mmtbx.validation_summary model.pdb

Prints a brief summary of validation criteria, including Ramachandran
statistics, rotamer outliers, clashscore, C-beta deviations, plus R-factors
and RMS(bonds)/RMS(angles) if found in PDB header.  (This is primarily used
for evaluating the output of refinement tests; general users are advised to
run phenix.model_vs_data or the validation GUI.)
""")
  parser = optparse.OptionParser()
  options, args = parser.parse_args(args)
  pdb_file = args[0]
  if (not os.path.isfile(pdb_file)):
    raise Sorry("Not a file: %s" % pdb_file)
  pdb_in = iotbx.pdb.input(pdb_file)
  hierarchy = pdb_in.construct_hierarchy()
  xrs = pdb_in.xray_structures_simple()
  crystal_symmetry=pdb_in.crystal_symmetry()
  if not crystal_symmetry:
    raise Sorry("Need crystal_symmetry in input PDB file")
  s = None
  extra = ""
  if (len(xrs) == 1):
    s = summary(pdb_file=pdb_file,crystal_symmetry=crystal_symmetry)
  else :
    s = ensemble(pdb_hierarchy=hierarchy,
      n_models=len(xrs),
      crystal_symmetry=crystal_symmetry)
    extra = " (%d models)" % len(xrs)
  print("", file=out)
  print("Validation summary for %s%s:" % (pdb_file, extra), file=out)
  s.show(out=out, prefix="  ", show_percentiles=True)
  print("", file=out)
  return s

if (__name__ == "__main__"):
  run(sys.argv[1:])
