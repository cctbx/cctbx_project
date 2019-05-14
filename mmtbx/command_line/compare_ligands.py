
from __future__ import absolute_import, division, print_function
from libtbx.utils import Usage, Sorry
import sys

master_phil = """
  pdb_file = None
    .type = path
    .multiple = True
  ligand_code = None
    .type = str
  max_distance_between_centers_of_mass = 8.0
    .type = float
  exclude_hydrogens = True
    .type = bool
  verbose = False
    .type = bool
"""

def run(args, out=sys.stdout):
  if (len(args) == 0) or ("--help" in args):
    raise Usage("""\
mmtbx.compare_ligands model1.pdb model2.pdb ligand_code=LIG

Given a pair of PDB files and a ligand ID, identifies equivalent copies of the
ligand in each model and calculates RMSDs between them.  Used for validating
ligand fitting tools in Phenix.""")

  import iotbx.phil
  class _cmdline(iotbx.phil.process_command_line_with_files):
    def process_other(self, arg):
      if (len(arg) <= 3) and (arg.isalnum()):
        return iotbx.phil.parse("""ligand_code=%s""" % arg)
  cmdline = _cmdline(
    args=args,
    master_phil_string=master_phil,
    pdb_file_def="pdb_file")
  params = cmdline.work.extract()
  if (len(params.pdb_file) != 2):
    raise Sorry("Exactly two PDB files required as input.")
  if (params.ligand_code is None):
    raise Sorry("Must specify ligand ID (ligand_code=LIG)")
  from mmtbx.validation import ligands
  rmsds, pbss = ligands.compare_ligands(
    ligand_code=params.ligand_code,
    pdb_file_1=params.pdb_file[0],
    pdb_file_2=params.pdb_file[1],
    exclude_hydrogens=params.exclude_hydrogens,
    verbose=params.verbose,
    out=out)
  if (len(rmsds) == 0):
    raise Sorry("No matching ligands found!")

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
