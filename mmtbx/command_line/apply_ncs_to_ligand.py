
from __future__ import absolute_import, division, print_function
from libtbx.str_utils import make_header
import sys

master_phil_str = """
ligand_code = LIG
  .type = str
atom_selection = None
  .type = atom_selection
add_to_model = True
  .type = bool
output_file = None
  .type = path
output_map = None
  .type = path
%s
"""

def run(args, out=sys.stdout):
  from mmtbx.command_line import load_model_and_data
  import mmtbx.ncs.ligands
  cmdline = load_model_and_data(
    args=args,
    master_phil=master_phil_str % mmtbx.ncs.ligands.ncs_ligand_phil,
    out=out,
    process_pdb_file=True,
    generate_input_phil=True,
    usage_string="""\
mmtbx.apply_ncs_to_ligand model.pdb data.mtz ligand_code=LIG ...

Given a multi-chain PDB file and a ligand residue name, find copies of the
ligand in the input file, identify NCS operators relating macromolecule chains,
and search for additional ligand sites by applying these operators.  Used to
complete ligand placement in cases where LigandFit (etc.) is only partially
successful.
""")
  pdb_hierarchy = cmdline.pdb_hierarchy
  fmodel = cmdline.fmodel
  params = cmdline.params
  if (params.output_file is None):
    params.output_file = "ncs_ligands.pdb"
  if (params.output_map is None):
    params.output_map = "ncs_ligands.mtz"
  make_header("Finding ligands by NCS operators", out=out)
  result = mmtbx.ncs.ligands.apply_ligand_ncs(
    pdb_hierarchy=pdb_hierarchy,
    fmodel=fmodel,
    params=params,
    ligand_code=params.ligand_code,
    atom_selection=None,
    add_new_ligands_to_pdb=params.add_to_model,
    log=out)
  result.write_pdb(params.output_file)
  result.write_maps(params.output_map)
  return result

if (__name__ == "__main__"):
  run(sys.argv[1:])
