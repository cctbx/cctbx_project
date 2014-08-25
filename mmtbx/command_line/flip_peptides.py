
from __future__ import division
import sys, os

#-----------------------------------------------------------------------
# XXX this isn't really very useful as a user-space program - it is
# mostly just here for testing purposes
def run (args, out=sys.stdout) :
  from libtbx.utils import Usage
  if (len(args) == 0) :
    print "Warning: this is not intended for general use."
    raise Usage("mmtbx.flip_peptides [model.pdb] [data.mtz] [params ...]")
  from mmtbx.refinement import flip_peptides, print_statistics
  import mmtbx.command_line
  master_phil = mmtbx.command_line.generate_master_phil_with_inputs(
    enable_automatic_twin_detection=True,
    phil_string="""
      flip_peptides {
        include scope mmtbx.refinement.flip_peptides.master_params_str
      }""")
  cmdline = mmtbx.command_line.load_model_and_data(
    args=args,
    master_phil=master_phil,
    out=out)
  print >> out, ""
  print_statistics.make_header("Analyzing peptide bonds", out=out)
  fmodel = cmdline.fmodel
  flip_peptides.run(
    fmodel=fmodel,
    geometry_restraints_manager=cmdline.geometry,
    pdb_hierarchy=cmdline.pdb_hierarchy,
    solvent_selection=None,
    params=cmdline.params.flip_peptides)
  if (len(cmdline.params.input.pdb.file_name) == 1) :
    pdb_file = cmdline.params.input.pdb.file_name[0]
    pdb_out = os.path.splitext(os.path.basename(pdb_file))[0] + "_new.pdb"
  else :
    pdb_out = "flipped.pdb"
  open(pdb_out, "w").write(cmdline.pdb_hierarchy.as_pdb_string())
  print >> out, ""
  print "wrote new PDB file to %s" % pdb_out

if __name__ == "__main__" :
  run(sys.argv[1:])
