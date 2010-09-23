
import sys, os

#-----------------------------------------------------------------------
# XXX this isn't really very useful as a user-space program - it is
# mostly just here for testing purposes
def run (args, out=sys.stdout) :
  from libtbx.utils import Sorry, Usage
  if (len(args) == 0) :
    print "Warning: this is not intended for general use."
    raise Usage("mmtbx.flip_peptides [model.pdb] [data.mtz] [params ...]")
  from mmtbx.refinement import flip_peptides, print_statistics
  from mmtbx import utils
  import iotbx.phil
  import libtbx.phil.command_line
  from cStringIO import StringIO
  master_phil = iotbx.phil.parse("""
    %s
    flip_peptides {
      include scope mmtbx.refinement.flip_peptides.master_params_str
    }
  """ % utils.cmdline_input_phil_str, process_includes=True)
  cmdline = utils.cmdline_load_pdb_and_data(
    args=args,
    master_phil=master_phil,
    out=out)
  print >> out, ""
  print_statistics.make_header("Analyzing peptide bonds", out=out)
  flip_peptides.run(
    fmodel=cmdline.fmodel,
    geometry_restraints_manager=cmdline.geometry,
    pdb_hierarchy=cmdline.pdb_hierarchy,
    solvent_selection=None)
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
