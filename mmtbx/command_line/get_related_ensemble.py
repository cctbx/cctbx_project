
from __future__ import absolute_import, division, print_function
from libtbx.str_utils import make_sub_header
from libtbx.utils import null_out
import sys

master_phil_str = """
model = None
  .type = path
chain_id = None
  .type = str
sequence = None
  .type = path
search_directory = None
  .type = path
output_file = ensemble.pdb
  .type = path
nproc = Auto
  .type = int
include scope mmtbx.building.make_library.master_phil
exclude_ids = None
  .type = strings
dry_run = False
  .type = bool
"""

def run(args, out=sys.stdout):
  from mmtbx.building import make_library
  import iotbx.phil
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil_string=master_phil_str,
    pdb_file_def="model",
    seq_file_def="sequence",
    directory_def="search_directory",
    usage_string="mmtbx.get_related_ensemble [model.pdb] [seq.fa] [...]")
  params = cmdline.work.extract()
  sequence = None
  if (params.model is None):
    raise Sorry("No model (PDB or mmCIF file) was specified.")
  if (params.sequence is not None):
    seq_file = cmdline.get_file(params.sequence, force_type="seq")
    n_seqs = len(seq_file.file_object)
    if (n_seqs > 1):
      print("%d sequences in file - will only use the first" % n_seqs, file=out)
    sequence = seq_file.file_object[0].sequence
  pdb_file = cmdline.get_file(params.model, force_type="pdb")
  hierarchy = pdb_file.file_object.hierarchy
  reference_hierarchy = iotbx.pdb.hierarchy.root()
  model = iotbx.pdb.hierarchy.model()
  reference_hierarchy.append_model(model)
  for chain in hierarchy.models()[0].chains():
    if (params.chain_id is None) or (chain.id == params.chain_id):
      if (not chain.is_protein()):
        if (chain.id == params.chain_id):
          print("warning: matching chain '%s' is not protein, skipping" % \
            chain.id, file=out)
        continue
      else :
        # TODO select based on sequence if provided
        new_chain = iotbx.pdb.hierarchy.chain(id=chain.id)
        model.append_chain(new_chain)
        # get rid of alternate conformations
        for residue_group in chain.residue_groups():
          atom_group = residue_group.atom_groups()[0]
          if (not atom_group.altloc.strip() in ['', 'A']):
            continue
          new_rg = iotbx.pdb.hierarchy.residue_group(
            resseq=residue_group.resseq,
            icode=residue_group.icode)
          new_ag = atom_group.detached_copy()
          new_ag.altloc = ''
          new_rg.append_atom_group(new_ag)
          new_chain.append_residue_group(new_rg)
        if (sequence is None):
          sequence = chain.as_padded_sequence(pad='X')
          print("Using sequence of chain '%s' (approx. %d residues)" % \
            (chain.id, len(sequence)), file=out)
        break
  if (sequence is None):
    raise Sorry("No protein sequence could be extracted based on these inputs.")
  make_sub_header("Finding related models and generating ensemble", out=out)
  ensemble = make_library.extract_and_superpose(
    reference_hierarchy=reference_hierarchy,
    search_directory=params.search_directory,
    sequence=sequence,
    params=params,
    out=out)
  f = null_out()
  if (params.output_file is not None):
    f = open(params.output_file, "w")
  print("Assembling moved models:", file=out)
  ensemble_hierarchy = ensemble.as_multi_model_hierarchy()
  for k in ensemble.selection_moved :
    source_info = ensemble.related_chains[k].source_info
    print("  Model %d: %s:%s" % (k+1, source_info,
      ensemble.related_chains[k].chain_id), file=out)
    f.write("REMARK model %d is from %s\n" % (k+1, source_info))
  f.write(ensemble_hierarchy.as_pdb_string())
  f.close()
  return ensemble_hierarchy

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
