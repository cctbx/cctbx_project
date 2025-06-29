"""Rearrange non-macromolecular heteroatom"""
# LIBTBX_SET_DISPATCHER_NAME phenix.sort_hetatms

# TODO sort waters by distance from macromolecule
# FIXME special treatment required for carbohydrates?

from __future__ import absolute_import, division, print_function
from libtbx.utils import Sorry, Usage, null_out
from libtbx import adopt_init_args
from libtbx import runtime_utils
import libtbx.phil
import operator
import os
import sys
from six.moves import zip

sorting_params_str = """
preserve_chain_id = False
  .type = bool
  .help = The default behavior is to group heteroatoms with the nearest \
    macromolecule chain, whose ID is inherited.  This parameter disables \
    the change of chain ID, and preserves the original chain ID.
waters_only = False
  .type = bool
  .help = Rearrange waters, but leave all other ligands alone.
sort_waters_by = none *b_iso
  .type = choice
  .help = Ordering of waters - by default it will sort them by the isotropic \
    B-factor.
  .caption = Don't_sort Isotropic_B
set_hetatm_record = True
  .type = bool
  .help = Convert ATOM to HETATM where appropriate.
  .short_caption = Set HETATM label
ignore_selection = None
  .type = atom_selection
  .help = Selection of atoms to skip.  Any residue group which overlaps with \
    this selection will be preserved with the original chain ID and numbering.
renumber = True
  .type = bool
  .help = Renumber heteroatoms once they are in new chains.
sequential_numbering = True
  .type = bool
  .help = If True, the heteroatoms will be renumbered starting from the next \
    available residue number after the end of the associated macromolecule \
    chain.  Otherwise, numbering will start from 1.
distance_cutoff = 6.0
  .type = float
  .help = Cutoff for identifying nearby macromolecule chains.  This should be \
    kept relatively small for speed reasons, but it may miss waters that are \
    far out in solvent channels.
  .input_size = 80
remove_waters_outside_radius = False
  .type = bool
  .help = Remove waters more than the specified distnace cutoff to the \
    nearest polymer chain (to avoid PDB complaints).
  .short_caption = Max. water-to-polymer distance
loose_chain_id = X
  .type = str
  .help = Chain ID assigned to heteroatoms that can't be mapped to a nearby \
    macromolecule chain.
  .short_caption = Loose chain ID
  .input_size = 48
model_skip_expand_with_mtrix = False
  .type = bool
  .help = If set to false then expand the model into NCS copies using \
    matrix (MTRIX) records present in PDB file
"""

master_params = """
file_name = None
  .type = path
  .short_caption = Model file
  .help = Input file
  .style = bold file_type:pdb input_file OnChange:extract_symmetry
output_file = None
  .type = path
  .style = bold file_type:pdb output_file new_file
unit_cell = None
  .type = unit_cell
space_group = None
  .type = space_group
ignore_symmetry = False
  .type = bool
  .help = Don't take symmetry-related chains into account when determining \
    the nearest macromolecule.
preserve_remarks = False
  .type = bool
  .help = Propagate all REMARK records to the output file.
verbose = False
  .type = bool
remove_hetatm_ter_records = True
  .type = bool
  .short_caption = Remove TER records after HETATM chains
  .help = The official PDB format only allows TER records at the end of \
    polymer chains, whereas the CCTBX PDB-handling tools will insert TER \
    after each chain of any type.  If this parameter is True, the extra \
    TER records will be removed.
%s
include scope libtbx.phil.interface.tracking_params
""" % sorting_params_str

def master_phil():
  return libtbx.phil.parse(master_params)

def sort_hetatms(
    pdb_hierarchy,
    xray_structure,
    params=None,
    verbose=False,
    return_pdb_hierarchy=True,
    log=null_out()):
  """
  Rearrange a PDB hierarchy so that heteroatoms are grouped with the closest
  macromolecule, accounting for symmetry.  Assorted ligands will be arranged
  first in each new chain, followed by waters.  See the PHIL block
  sorting_params_str for options.
  """
  if (params is None):
    import iotbx.phil
    params = iotbx.phil.parse(sorting_params_str).fetch().extract()
  from iotbx import pdb
  from scitbx.array_family import flex
  pdb_atoms = pdb_hierarchy.atoms()
  pdb_atoms.reset_i_seq()
  # the i_seq will be reset in various places, so I substitute the tmp
  # attribute
  for atom in pdb_atoms :
    atom.tmp = atom.i_seq # XXX why can't this be an array operation?
  pdb_atoms = pdb_hierarchy.deep_copy().atoms()
  n_atoms = pdb_atoms.size()
  assert (n_atoms == len(xray_structure.scatterers()))
  ignore_selection = flex.bool(n_atoms, False)
  mm_selection = flex.bool(n_atoms, False)
  hd_selection = xray_structure.hd_selection()
  sites_frac = xray_structure.sites_frac()
  new_sites_cart = xray_structure.sites_cart()
  pair_asu_table = xray_structure.pair_asu_table(
    distance_cutoff=params.distance_cutoff)
  asu_mappings = pair_asu_table.asu_mappings()
  asu_table = pair_asu_table.table()
  mm_chains = []
  hetatm_chains = []
  if (params.ignore_selection is not None):
    sel_cache = pdb_hierarchy.atom_selection_cache()
    ignore_selection = sel_cache.selection(params.ignore_selection)
  assert (len(pdb_hierarchy.models()) == 1)
  for chain in pdb_hierarchy.only_model().chains():
    chain_atoms = chain.atoms()
    het_sel = chain_atoms.extract_hetero()
    if (het_sel.size() == chain_atoms.size()):
      hetatm_chains.append(chain)
    elif ((chain.is_protein(ignore_water=False)) or
          (chain.is_na(ignore_water=False))):
      mm_chains.append(chain)
      i_seqs = chain_atoms.extract_tmp_as_size_t()
      mm_selection.set_selected(i_seqs, True)
    else :
      hetatm_chains.append(chain)
  if (len(hetatm_chains) == 0):
    print("No heteroatoms - hierarchy will not be modified.", file=log)
    if (return_pdb_hierarchy):
      return pdb_hierarchy
    else :
      return sort_hetatms_result(
        pdb_hierarchy=pdb_hierarchy,
        n_mm_chains=len(mm_chains),
        n_het_residues=0)
  n_het_residues = 0
  new_hierarchy = pdb.hierarchy.root()
  new_model = pdb.hierarchy.model()
  new_chain_i_seqs = []
  new_hierarchy.append_model(new_model)
  new_hetatm_chains = []
  new_hetatm_chain_ids = []
  new_start_resseq = []
  for chain in mm_chains :
    new_chain_i_seqs.append(flex.size_t(chain.atoms().extract_tmp_as_size_t()))
    new_model.append_chain(chain.detached_copy())
    new_hetatm_chains.append(pdb.hierarchy.chain(id=chain.id))
    if (params.preserve_chain_id) and (chain.id in new_hetatm_chain_ids):
      print("Warning: chain ID '%s' is duplicated", file=log)
    new_hetatm_chain_ids.append(chain.id)
    if (params.sequential_numbering):
      last_resseq = chain.residue_groups()[-1].resseq_as_int()
      new_start_resseq.append(last_resseq + 1)
    else :
      new_start_resseq.append(1)
  hetatm_residue_groups = []
  hetatm_residue_chain_ids = []
  loose_chain_id = params.loose_chain_id
  if (loose_chain_id is None):
    loose_chain_id = " "
  loose_residues = pdb.hierarchy.chain(id=loose_chain_id)
  preserve_chains = []
  for chain in hetatm_chains :
    chain_id = chain.id
    for rg in chain.residue_groups():
      n_het_residues += 1
      if (params.preserve_chain_id):
        if (chain_id in new_hetatm_chain_ids):
          i_chain = new_hetatm_chain_ids.index(chain_id)
          new_hetatm_chains[i_chain].append_residue_group(rg.detached_copy())
        else :
          print("Warning: no corresponding macromolecule chain match for %s %s" % \
            (chain_id, rg.resid()), file=log)
          loose_residues.append_residue_group(rg.detached_copy())
        chain.remove_residue_group(rg)
        continue
      keep_residue_group = False
      rg_atoms = rg.atoms()
      i_seqs = rg_atoms.extract_tmp_as_size_t()
      atom_groups = rg.atom_groups()
      if ((params.waters_only) and
          (not atom_groups[0].resname in ["HOH","WAT"])):
        keep_residue_group = True
      else :
        for i_seq in i_seqs :
          if (ignore_selection[i_seq]):
            keep_residue_group = True
            break
      if (not keep_residue_group):
        rg_copy = rg.detached_copy()
        for new_atom, old_atom in zip(rg_copy.atoms(), rg_atoms):
          new_atom.tmp = old_atom.tmp # detached_copy() doesn't preserve tmp
        hetatm_residue_groups.append(rg_copy)
        hetatm_residue_chain_ids.append(chain_id)
        chain.remove_residue_group(rg)
    if (len(chain.residue_groups()) > 0):
      preserve_chains.append(chain.detached_copy())
  for chain in preserve_chains :
    new_model.append_chain(chain)
  unit_cell = xray_structure.unit_cell()
  n_deleted = 0
  if (not params.preserve_chain_id):
    for k, rg in enumerate(hetatm_residue_groups):
      chain_id = hetatm_residue_chain_ids[k]
      atom_groups = rg.atom_groups()
      if (len(atom_groups) == 0):
        continue
      is_water = (atom_groups[0].resname in ["HOH", "WAT", "DOD"])
      rg_atoms = rg.atoms()
      i_seqs = rg_atoms.extract_tmp_as_size_t()
      closest_distance = sys.maxsize
      closest_i_seq = None
      closest_rt_mx = None
      for i_seq, atom in zip(i_seqs, rg_atoms):
        if (params.set_hetatm_record):
          atom.hetero = True
        if (hd_selection[i_seq]):
          continue
        site_i = sites_frac[i_seq]
        asu_dict = asu_table[i_seq]
        rt_mx_i_inv = asu_mappings.get_rt_mx(i_seq, 0).inverse()
        for j_seq, j_sym_groups in asu_dict.items():
          if (hd_selection[j_seq]) or (not mm_selection[j_seq]):
            continue
          site_j = sites_frac[j_seq]
          for j_sym_group in j_sym_groups:
            rt_mx = rt_mx_i_inv.multiply(asu_mappings.get_rt_mx(j_seq,
              j_sym_group[0]))
            site_ji = rt_mx * site_j
            dxyz = unit_cell.distance(site_i, site_ji)
            if (dxyz < closest_distance):
              closest_distance = dxyz
              closest_rt_mx = rt_mx.inverse() # XXX I hope this is right...
              closest_i_seq = j_seq
      if (closest_i_seq is None):
        # XXX possible bug: what about waters H-bonded to ligands, but still
        # outside radius of macromolecule?
        if (is_water and params.remove_waters_outside_radius):
          print("Water %s is not near any polymer chain, will delete" \
            % rg.id_str(), file=log)
          n_deleted += len(rg.atoms())
          continue
        # XXX not sure what to do here - if we have only one macromolecule
        # chain, does it make more sense to keep all hetatms grouped together
        # regardless of distance?  e.g. waters in 1BS2
        else :
          print("Residue group %s %s is not near any macromolecule chain" % \
            (chain_id, rg.resid()), file=log)
          loose_residues.append_residue_group(rg)
      else :
        for j_seqs, hetatm_chain in zip(new_chain_i_seqs, new_hetatm_chains):
          if (closest_i_seq in j_seqs):
            if (verbose):
              if (closest_rt_mx.is_unit_mx()):
                print("Residue group %s added to chain %s (distance = %.3f)" % \
                  (rg.atoms()[0].id_str(), hetatm_chain.id, closest_distance), file=log)
              else :
                print(("Residue group %s added to chain %s "+
                   "(distance = %.3f, symop = %s)") % \
                  (rg.atoms()[0].id_str(), hetatm_chain.id, closest_distance,
                   str(closest_rt_mx)), file=log)
            if (not closest_rt_mx.is_unit_mx()):
              # closest macromolecule is in another ASU, so map the hetatms to
              # be near the copy in the current ASU
              for atom in rg.atoms():
                site_frac = unit_cell.fractionalize(site_cart=atom.xyz)
                new_site_frac = closest_rt_mx * site_frac
                atom.xyz = unit_cell.orthogonalize(site_frac=new_site_frac)
            hetatm_chain.append_residue_group(rg)
            break
        else :
          raise RuntimeError("Can't find chain for i_seq=%d" % closest_i_seq)
  # even if waters aren't sorted, we still want them to come last
  for chain in new_hetatm_chains :
    waters_and_b_iso = []
    for rg in chain.residue_groups():
      ags = rg.atom_groups()
      if (ags[0].resname in ["WAT","HOH"]):
        b_iso = flex.mean(rg.atoms().extract_b())
        waters_and_b_iso.append((rg, b_iso))
        chain.remove_residue_group(rg)
    if (len(waters_and_b_iso) > 0):
      if (params.sort_waters_by != "none"):
        waters_and_b_iso.sort(key=operator.itemgetter(1))
      for water, b_iso in waters_and_b_iso :
        chain.append_residue_group(water)
  if (params.renumber):
    for chain, start_resseq in zip(new_hetatm_chains, new_start_resseq):
      resseq = start_resseq
      for rg in chain.residue_groups():
        rg.resseq = resseq
        resseq += 1
  for chain in new_hetatm_chains :
    for residue_group in chain.residue_groups():
      residue_group.link_to_previous = True # suppress BREAK records
    if len(chain.atoms()) > 0:  # Skip empty chains
      new_model.append_chain(chain)
  if (len(loose_residues.residue_groups()) > 0):
    new_model.append_chain(loose_residues)
  n_atoms_new = len(new_hierarchy.atoms())
  if (n_atoms_new != n_atoms - n_deleted):
    raise RuntimeError("Atom counts do not match: %d --> %d" % (n_atoms,
      n_atoms_new))
  if (n_deleted > 0):
    print("WARNING: %d atoms removed from model" % n_deleted, file=log)
    print("  You must refine this model again before deposition!", file=log)
  if (return_pdb_hierarchy):
    return new_hierarchy
  else :
    return sort_hetatms_result(
      pdb_hierarchy=new_hierarchy,
      n_mm_chains=len(mm_chains),
      n_het_residues=n_het_residues)

def run(args, out=sys.stdout, sorting_params=None):
  import iotbx.phil
  if (len(args) == 0) or ("--help" in args):
    raise Usage("""\
mmtbx.sort_hetatms model.pdb [options]

Rearrange non-macromolecular heteroatoms (waters and other ligands, by default)
into new chains having the same ID as the closest macromolecule chain.  Will
also renumber residues and sort waters by B-factor with default settings.

Full parameters:

%s
""" % iotbx.phil.parse(master_params).as_str(prefix="  "))
  import iotbx.pdb
  from cctbx import crystal
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil_string=master_params,
    pdb_file_def="file_name")
  params = cmdline.work.extract()
  validate_params(params)
  pdb_in = iotbx.pdb.input(params.file_name)
  pdb_symm = pdb_in.crystal_symmetry()
  space_group = params.space_group
  unit_cell = params.unit_cell
  if (pdb_symm is None) and (not params.ignore_symmetry):
    if (space_group is None) or (unit_cell is None):
      raise Sorry("Crystal symmetry information is required; please specify "+
        "the space_group and unit_cell parameters.")
  else :
    if (space_group is None) and (pdb_symm is not None):
      space_group = pdb_symm.space_group_info()
    if (unit_cell is None) and (pdb_symm is not None):
      unit_cell = pdb_symm.unit_cell()
  final_symm = None
  if (not params.ignore_symmetry):
    final_symm = crystal.symmetry(
      space_group_info=space_group,
      unit_cell=unit_cell)
  pdb_hierarchy = pdb_in.construct_hierarchy()
  xray_structure = pdb_in.xray_structure_simple(
    crystal_symmetry=final_symm)
  if (sorting_params is None):
    sorting_params = params
  result = sort_hetatms(
    pdb_hierarchy=pdb_hierarchy,
    xray_structure=xray_structure,
    params=sorting_params,
    verbose=params.verbose,
    return_pdb_hierarchy=False,
    log=out)
  if (params.output_file is None):
    params.output_file = os.path.splitext(
      os.path.basename(params.file_name))[0] + "_sorted.pdb"
  if (not result.pdb_hierarchy.fits_in_pdb_format()):
    params.output_file = result.pdb_hierarchy.write_pdb_or_mmcif_file(
      target_filename = params.output_file)

  else: # standard pdb file
    f = open(params.output_file, "w")
    if (params.preserve_remarks):
      remarks = pdb_in.remark_section()
      if (len(remarks) > 0):
        f.write("\n".join(remarks))
        f.write("\n")
    pdb_str = result.pdb_hierarchy.as_pdb_string(crystal_symmetry=final_symm)
    if (params.remove_hetatm_ter_records):
      n_hetatm = n_atom = 0
      for line in pdb_str.splitlines():
        if (line[0:3] == "TER"):
          if (n_atom != 0):
            f.write("%s\n" % line)
          n_atom = n_hetatm = 0
          continue
        elif (line.startswith("HETATM")):
          n_hetatm += 1
        elif (line.startswith("ATOM")):
          n_atom += 1
        f.write("%s\n" % line)
    else :
      f.write(pdb_str)
    f.write("END")
    f.close()

  print("Wrote %s" % params.output_file, file=out)
  out.flush()
  return sort_hetatms_result(
    file_name=os.path.abspath(params.output_file),
    n_mm_chains=result.n_mm_chains,
    n_het_residues=result.n_het_residues)

class sort_hetatms_result(object):
  def __init__(self, n_mm_chains, n_het_residues, pdb_hierarchy=None,
      file_name=None):
    adopt_init_args(self, locals())
    assert ([pdb_hierarchy, file_name].count(None) == 1)

  def finish_job(self):
    return ([(self.file_name, "Modified model")],
            [("Number of macromolecule chains", self.n_mm_chains),
             ("Number of heteroatom residues", self.n_het_residues)])

def validate_params(params):
  if (params.file_name is None):
    raise Sorry("Model file (file_name) not specified.")
  if not os.path.isfile(params.file_name):
    raise Sorry("Model file (%s) is missing." %(params.file_name))
  from iotbx.data_manager import DataManager
  if (params.model_skip_expand_with_mtrix):
    # phenix.famos doesn't need NCS expanded model and some files in PDB contain
    # faulty matrices causing a Sorry on expansion.
    # Use this keyword for phenix.famos to avoid expansion
    dm = DataManager(custom_options=['model_skip_expand_with_mtrix'])
  else:
    dm = DataManager()
  m = dm.get_model(params.file_name)
  if not m:
    raise Sorry("Unable to read the file %s" %(params.file_name))

  return True

class launcher(runtime_utils.target_with_save_result):
  def run(self):
    return run(args=list(self.args), out=sys.stdout)


if (__name__ == "__main__"):
  run(sys.argv[1:])

