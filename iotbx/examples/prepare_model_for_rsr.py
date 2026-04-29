"""Preparation of model for real-space refinement"""
from __future__ import absolute_import, division, print_function
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out, Sorry
import sys

"""
Example: preparation of model for real-space refinement:
  - Expand MTRIX
  - Expand BIOMT
  - Build geometry restraints manager
  - Show initial statistics.
Should support PDB and mmCIF formats.

Usage:
  python prepare_model_for_rsr.py <model>
"""

def show_ss_counts(model):
  ss = model.get_ss_annotation()
  print("Number of helices:", ss.get_n_helices())
  print("Number of sheets :", ss.get_n_sheets())



def run(args):
  assert len(args) == 1
  # Read file into pdb_input class
  inp = iotbx.pdb.input(file_name=args[0])

  # create a model manager
  # Catch Sorry about MTRIX here.
  model = mmtbx.model.manager(
      model_input = inp,
      restraint_objects = None, # these are ligands if any [('fname', cif_object), ()]
      log = null_out(),
      )
  print("="*80)
  print("number of atoms with MTRIX multiplication:", model.get_number_of_atoms())
  show_ss_counts(model)

  # Expand with BIOMT if needed. MTRIX are already expanded by default
  # Catch case when both MTRIX and BIOMT present, or other Sorry raised by
  # BIOMT handling.
  # LIMITATION: this should be done before any selections made on model.manager
  double_counter = 0
  try:
    model.expand_with_BIOMT_records()
  except Sorry as e:
    if str(e).startswith("Model has been already expanded"):
      double_counter += 1
  print("="*80)
  print("number of atoms with BIOMT multiplication:", model.get_number_of_atoms())
  show_ss_counts(model)

  # Get default params
  pdb_int_params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  # Set whatever you want
  pdb_int_params.pdb_interpretation.secondary_structure.protein.enabled = True
  pdb_int_params.pdb_interpretation.ncs_search.enabled = True
  pdb_int_params.pdb_interpretation.ncs_search.residue_match_radius = 999
  pdb_int_params.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None

  #pdb_int_params.pdb_interpretation.nonbonded_weight = None

  # set the params. Note, that GRM would be dropped, even if it was already
  # constructed. In this example it is not yet constructed.
  model.set_pdb_interpretation_params(params=pdb_int_params)
  grm = model.get_restraints_manager()

  # Not clear which one should be used at the moment
  gs = model.geometry_statistics()
  gs.show()
  # The second way
  msi = model.get_model_statistics_info()
  msi.show_remark_3()

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
