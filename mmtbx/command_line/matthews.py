
# TODO tests

from __future__ import absolute_import, division, print_function
from iotbx import crystal_symmetry_from_any
import iotbx.bioinformatics
import iotbx.phil
from cctbx import crystal
from libtbx.utils import Sorry
import sys

master_phil_str = """
data = None
  .type = path
labels = None
  .type = str
sequence = None
  .type = path
model = None
  .type = path
space_group = None
  .type = space_group
unit_cell = None
  .type = unit_cell
n_residues = None
  .type = int(value_min=1)
  .optional = True
n_bases = None
  .type = int(value_min=1)
  .optional = True
"""

def run(args, out=sys.stdout):
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil_string=master_phil_str,
    pdb_file_def="model",
    reflection_file_def="data",
    seq_file_def="sequence",
    space_group_def="space_group",
    unit_cell_def="unit_cell",
    integer_def="n_residues",
    usage_string="""\
phenix.matthews [data.hkl] [space_group] [unit_cel] [sequence] [n_residues] ...

Calculate the expected Matthews coefficient given the crystal symmetry and
crystallized molecule(s).
""")
  params = cmdline.work.extract()
  if (params.space_group is None) or (params.unit_cell is None):
    if (params.data is None):
      raise Sorry("You must supply both a space group and a unit cell (or "+
        "a data file containing this information).")
    else :
      symm = crystal_symmetry_from_any.extract_from(file_name=params.data)
      space_group_from_file = symm.space_group()
      if (params.space_group is None):
        if (space_group_from_file is not None):
          params.space_group = symm.space_group()
      elif (space_group_from_file is not None):
        if (space_group_from_file != params.space_group):
          print("WARNING: space group mismatch between command line "+\
            "and file:", file=out)
          print("  %s (cmdline), %s (file)" % (params.space_group,
            space_group_from_file), file=out)
      if (params.unit_cell is None):
        params.unit_cell = symm.unit_cell()
  validate_params(params, check_symmetry=True)
  if (params.sequence is not None):
    assert (params.n_residues == params.n_bases == None)
    seq_comp = iotbx.bioinformatics.composition_from_sequence_file(
      file_name=params.sequence,
      log=out)
    if (seq_comp is not None):
      params.n_residues = seq_comp.n_residues
      params.n_bases = seq_comp.n_bases
    else :
      raise Sorry("No composition information could be obtained from the "+
        "sequence file.")
  elif (params.model is not None):
    assert (params.n_residues == params.n_bases == None)
    from iotbx.file_reader import any_file
    params.n_residues = 0
    params.n_bases = 0
    pdb_in = any_file(params.model)
    hierarchy = pdb_in.file_object.hierarchy
    for chain in hierarchy.models()[0].chains():
      if chain.is_protein():
        params.n_residues += chain.residue_groups_size()
      elif chain.is_na():
        params.n_bases += chain.residue_groups_size()
  print("Space group: %s" % params.space_group, file=out)
  print("Unit cell: %s" % params.unit_cell, file=out)
  if (params.n_residues > 0):
    print("Number of residues: %d" % params.n_residues, file=out)
  if (params.n_bases > 0):
    print("Number of bases: %d" % params.n_bases, file=out)
  symm = crystal.symmetry(
    space_group_info=params.space_group,
    unit_cell=params.unit_cell)
  from mmtbx.scaling import matthews
  result = matthews.matthews_rupp(
    crystal_symmetry=symm,
    n_residues=params.n_residues,
    n_bases=params.n_bases)
  result.show(out=out)
  return result

def validate_params(params, check_symmetry=False):
  if (params.sequence is None) and (params.model is None):
    if (params.n_residues is None) and (params.n_bases is None):
      pass
      #raise Sorry("You must specify the composition of the crystallized "+
      #  "entity - either a sequence, a partial model, or the number of "+
      #  "protein residues or nucleic acid bases.")
  else :
    if (params.n_residues is not None) or (params.n_bases is not None):
      raise Sorry("You may only specify a sequence file OR a partial model "+
        "OR the numbers of residues and/or bases.")
  if (check_symmetry):
    if (params.space_group is None) or (params.unit_cell is None):
      raise Sorry("You must supply both a space group and a unit cell (or "+
        "a data file containing this information).")
  return True

if (__name__ == "__main__"):
  run(sys.argv[1:])
