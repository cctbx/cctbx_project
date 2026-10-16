"""
probe2 output.cxx_dots and clashscore2's per-atom-pair contact records must give the same
results as the Python dot loop and the per-dot JSON parse they replace.
"""
from __future__ import absolute_import, division, print_function
import re, time
import iotbx.cli_parser
import iotbx.phil
import iotbx.pdb
import libtbx.load_env
import mmtbx.model
from iotbx.data_manager import DataManager
from libtbx.utils import null_out
from mmtbx.hydrogens import reduce_hydrogen
from mmtbx.probe import Helpers
from mmtbx.programs import probe2
from mmtbx.validation import clashscore2

# Alternates on a backbone and a side chain, partial-occupancy waters (from tst_clashscore_2).
pdb_str = """\
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1
ATOM      1  N   ASN A   1      12.388  -5.203  29.298  1.00 13.23           N
ATOM      2  CA  ASN A   1      12.398  -3.931  28.500  1.00  6.85           C
ATOM      3  CB  ASN A   1      13.493  -3.033  29.043  1.00  7.32           C
ATOM      4  CG  ASN A   1      13.414  -1.559  28.658  1.00  7.93           C
ATOM      5  OD1 ASN A   1      14.388  -0.782  28.842  1.00 10.99           O
ATOM      6  ND2 ASN A   1      12.315  -1.201  28.134  1.00 11.12           N
ATOM      7  C  AASN A   1      12.384  -4.463  27.091  0.45  6.81           C
ATOM      8  O  AASN A   1      12.645  -5.635  26.750  0.45  8.28           O
ATOM      9  C  BASN A   1      12.891  -4.205  27.064  0.55  9.54           C
ATOM     10  O  BASN A   1      14.004  -4.638  26.750  0.55 12.04           O
ATOM    180  N   GLU A  25      17.119   0.317  19.886  1.00  4.46           N
ATOM    181  CA  GLU A  25      17.951   0.073  21.066  1.00  4.05           C
ATOM    182  C   GLU A  25      17.083  -0.029  22.303  1.00  3.92           C
ATOM    183  O   GLU A  25      16.053   0.642  22.379  1.00  4.03           O
ATOM    184  CB  GLU A  25      18.993   1.197  21.334  1.00  5.58           C
ATOM    185  CG  GLU A  25      20.050   1.356  20.259  1.00  6.70           C
ATOM    186  CD AGLU A  25      21.014   2.440  20.770  0.47  7.87           C
ATOM    187  OE1AGLU A  25      22.134   2.188  21.301  0.47  9.47           O
ATOM    188  OE2AGLU A  25      20.526   3.611  20.664  0.47  7.81           O
ATOM    189  CD BGLU A  25      20.877   2.608  20.405  0.53  8.00           C
ATOM    190  OE1BGLU A  25      21.950   2.605  19.769  0.53  8.34           O
ATOM    191  OE2BGLU A  25      20.372   3.528  21.137  0.53  7.48           O
ATOM    192  N   VAL A  26      17.559  -0.774  23.291  1.00  4.51           N
ATOM    193  CA  VAL A  26      16.931  -0.812  24.603  1.00  4.35           C
ATOM    194  C   VAL A  26      17.394   0.427  25.384  1.00  4.53           C
ATOM    195  O   VAL A  26      18.565   0.603  25.717  1.00  6.63           O
ATOM    196  CB  VAL A  26      17.285  -2.120  25.322  1.00  6.28           C
ATOM    197  CG1 VAL A  26      16.710  -2.090  26.737  1.00  8.38           C
ATOM    198  CG2 VAL A  26      16.793  -3.328  24.541  1.00  7.99           C
ATOM    199  N   VAL A  27      16.430   1.306  25.641  1.00  4.09           N
ATOM    200  CA  VAL A  27      16.630   2.593  26.282  1.00  4.04           C
ATOM    201  C   VAL A  27      15.546   2.821  27.325  1.00  3.81           C
ATOM    202  O   VAL A  27      14.515   2.145  27.331  1.00  4.70           O
ATOM    203  CB  VAL A  27      16.645   3.734  25.238  1.00  4.32           C
ATOM    204  CG1 VAL A  27      17.867   3.681  24.325  1.00  5.54           C
ATOM    205  CG2 VAL A  27      15.361   3.732  24.400  1.00  4.82           C
ATOM    206  N   THR A  28      15.797   3.825  28.171  1.00  3.67           N
ATOM    207  CA  THR A  28      14.793   4.268  29.126  1.00  3.78           C
ATOM    208  C   THR A  28      14.774   5.796  29.120  1.00  3.63           C
ATOM    209  O   THR A  28      15.846   6.424  29.117  1.00  4.10           O
ATOM    210  CB  THR A  28      15.041   3.745  30.540  1.00  4.19           C
ATOM    211  OG1 THR A  28      16.342   4.127  31.000  1.00  4.27           O
ATOM    212  CG2 THR A  28      14.895   2.227  30.627  1.00  5.21           C
ATOM    213  N   PRO A  29      13.578   6.399  29.122  1.00  3.53           N
ATOM    214  CA  PRO A  29      12.262   5.767  29.064  1.00  3.65           C
ATOM    215  C   PRO A  29      11.992   5.133  27.695  1.00  3.56           C
ATOM    216  O   PRO A  29      12.719   5.286  26.720  1.00  5.32           O
ATOM    217  CB  PRO A  29      11.293   6.907  29.372  1.00  4.74           C
ATOM    218  CG  PRO A  29      12.017   8.133  28.814  1.00  4.83           C
ATOM    219  CD  PRO A  29      13.477   7.862  29.201  1.00  4.02           C
ATOM    220  N   MET A  30      10.875   4.415  27.628  1.00  4.82           N
ATOM    221  CA  MET A  30      10.472   3.706  26.429  1.00  3.92           C
ATOM    222  C   MET A  30       9.962   4.670  25.386  1.00  3.81           C
ATOM    223  O   MET A  30       9.205   5.595  25.690  1.00  4.84           O
ATOM    224  CB  MET A  30       9.331   2.734  26.807  1.00  4.57           C
ATOM    225  CG  MET A  30       8.873   1.861  25.654  1.00  4.88           C
ATOM    226  SD  MET A  30      10.130   0.613  25.308  1.00  5.54           S
ATOM    227  CE AMET A  30       9.590   0.158  23.659  0.80  5.92           C
ATOM    228  CE BMET A  30       9.842  -0.500  26.625  0.20 14.35           C
ATOM    229  N   GLY A  31      10.375   4.449  24.142  1.00  3.74           N
ATOM    230  CA  GLY A  31       9.846   5.184  23.024  1.00  4.22           C
ATOM    231  C   GLY A  31       8.887   4.357  22.184  1.00  3.46           C
ATOM    232  O   GLY A  31       8.350   3.328  22.642  1.00  3.89           O
ATOM    233  N   ILE A  32       8.655   4.794  20.950  1.00  3.51           N
ATOM    234  CA  ILE A  32       7.877   3.976  20.038  1.00  3.38           C
ATOM    235  C   ILE A  32       8.559   2.602  19.944  1.00  3.44           C
ATOM    236  O   ILE A  32       9.782   2.531  19.754  1.00  3.63           O
ATOM    237  CB  ILE A  32       7.807   4.651  18.658  1.00  3.57           C
ATOM    238  CG1 ILE A  32       7.052   5.977  18.750  1.00  4.18           C
ATOM    239  CG2 ILE A  32       7.139   3.697  17.669  1.00  4.59           C
ATOM    240  CD1 ILE A  32       7.076   6.794  17.464  1.00  4.77           C
ATOM    241  N   PRO A  33       7.817   1.492  20.034  1.00  3.75           N
ATOM    242  CA  PRO A  33       8.459   0.197  19.927  1.00  4.43           C
ATOM    243  C   PRO A  33       9.162   0.012  18.587  1.00  3.94           C
ATOM    244  O   PRO A  33       8.668   0.406  17.535  1.00  4.50           O
ATOM    245  CB  PRO A  33       7.289  -0.810  20.093  1.00  6.13           C
ATOM    246  CG  PRO A  33       6.309  -0.027  20.992  1.00  6.43           C
ATOM    247  CD  PRO A  33       6.394   1.383  20.414  1.00  4.84           C
ATOM    248  N   ALA A  34      10.310  -0.693  18.648  1.00  4.36           N
ATOM    249  CA  ALA A  34      11.123  -0.882  17.448  1.00  4.68           C
ATOM    250  C   ALA A  34      10.374  -1.543  16.299  1.00  4.72           C
ATOM    251  O   ALA A  34      10.729  -1.285  15.141  1.00  5.19           O
ATOM    252  CB  ALA A  34      12.370  -1.700  17.814  1.00  7.22           C
TER
HETATM  156  O   HOH A 173       9.659  -2.091  27.019  0.29  6.50           O
HETATM  191  O   HOH A 208      24.121   2.918  18.439  0.20  7.03           O
HETATM  217  O   HOH A 234      10.755   0.622  29.147  0.23  8.84           O
HETATM  230  O   HOH A 247      23.962   1.544  18.100  0.26  4.79           O
END
"""

def hydrogenated():
  """(hierarchy, crystal symmetry) of the model with H placed."""
  inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  model = mmtbx.model.manager(model_input=inp, log=null_out())
  ph = reduce_hydrogen.place_hydrogens(model=model)
  ph.run()
  return ph.get_model().get_hierarchy(), model.crystal_symmetry()

def fresh_model(hierarchy, symmetry):
  # A new model each run: probe2 adds Phantom Hydrogens to the waters of the one it runs on.
  return mmtbx.model.manager(model_input=None, pdb_hierarchy=hierarchy.deep_copy(),
                             crystal_symmetry=symmetry, log=null_out())

def run_probe2(model, args):
  parser = iotbx.cli_parser.CCTBXParser(program_class=probe2.Program, logger=null_out())
  parser.parse_args(args + ["output.write_files=False", "ignore_lack_of_explicit_hydrogens=True"])
  dm = DataManager()
  dm.add_model("model", model)
  p2 = probe2.Program(dm, parser.working_phil.extract(), master_phil=parser.master_phil,
                      logger=null_out())
  p2.overrideModel(model)
  # The header carries the time of the run.
  return re.sub(r"run \d{4}-\d\d-\d\d \d\d:\d\d:\d\d", "run <time>", p2.run()[1])

def two_models(hierarchy):
  h = hierarchy.deep_copy()
  h.only_model().id = "1"
  second = h.only_model().detached_copy()
  second.id = "2"
  h.append_model(second)
  return h

RAW, KIN = ["output.format=raw"], ["output.format=kinemage"]

# One case per ContactDB input (set_params, atom flags, bonds) and per kinemage writer option,
# since kinemage without count_dots or contact_summary is written from the C++ arrays.
CXX_DOTS_CASES = [
  KIN, RAW, ["output.format=json"], ["output.count_dots=True"],
  # Waters as source take the Phantom Hydrogen path.
  ["approach=both", "source_selection=water", "target_selection=all"] + RAW,
  ["approach=both", "source_selection=protein", "target_selection=water"] + KIN,
  ["approach=once", "source_selection=protein", "target_selection=not protein"] + RAW,
  ["output.contact_summary=True"] + KIN,
  ["output.condensed=True"] + RAW,
  ["minimum_occupancy=0.5"] + RAW,
  ["excluded_bond_chain_length=3"] + RAW,
  ["excluded_bond_chain_length=5"] + RAW,
  ["include_mainchain_mainchain=False"] + RAW,
  ["include_water_water=True", "source_selection=all"] + RAW,
  ["overlap_scale_factor=0.3"] + RAW,
  ["probe.probe_radius=0.0"] + RAW,
  # Weak H-bond and worse overlap dots need these; the lower cutoffs give this model bad overlaps.
  ["probe.allow_weak_hydrogen_bonds=True", "output.separate_worse_clashes=True",
   "probe.clash_cutoff=0.1", "probe.worse_clash_cutoff=0.2"] + RAW,
  ["output.separate_worse_clashes=True"] + KIN,
  ["output.report_hydrogen_bonds=False", "output.report_clashes=False"] + RAW,
  ["output.report_vdws=False"] + KIN,
  ["output.bin_gaps=True"] + KIN,
  ["output.color_by_gap=False", "output.merge_contacts=False"] + KIN,
  ["output.atoms_are_masters=True", "output.altid_as_pointmaster=False"] + KIN,
]

def exercise_cxx_dots():
  hierarchy, symmetry = hydrogenated()
  # Two identical models: dots between them would be all clashes.
  cases = [(hierarchy, args) for args in CXX_DOTS_CASES] + [(two_models(hierarchy), RAW)]
  for h, args in cases:
    python = run_probe2(fresh_model(h, symmetry), args + ["output.cxx_dots=False"])
    cxx = run_probe2(fresh_model(h, symmetry), args + ["output.cxx_dots=True"])
    assert len(python.splitlines()) > 5, (args, python)
    assert cxx == python, args

def run_clashscore2(**kwds):
  dm = DataManager()
  dm.process_model_str("model", pdb_str)
  probe_params = iotbx.phil.parse(Helpers.probe_phil_parameters).extract()
  return clashscore2.clashscore2(probe_parameters=probe_params.probe, data_manager=dm,
                                 keep_hydrogens=False, out=null_out(), **kwds)

def exercise_contact_records():
  manager = clashscore2.probe_clashscore_manager
  assert manager.use_contact_records
  records = run_clashscore2()
  manager.use_contact_records = False
  try:
    dots = run_clashscore2()
  finally:
    manager.use_contact_records = True
  assert len(records.results) > 0
  assert records.as_JSON() == dots.as_JSON()
  assert records.get_clashscore() == dots.get_clashscore()

if (__name__ == "__main__"):
  if (libtbx.env.has_module("reduce") and libtbx.env.has_module("probe")):
    t0 = time.time()
    exercise_cxx_dots()
    exercise_contact_records()
    print("OK. Time: %8.3f" % (time.time() - t0))
  else:
    print("Skipping tst_probe2_contact_records: probe or reduce not configured")
