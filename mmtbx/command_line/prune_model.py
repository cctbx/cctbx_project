
# TODO trim sidechains one atom at a time

import libtbx.phil
from libtbx.str_utils import make_header
from libtbx.utils import Usage, multi_out
from libtbx import group_args
from cStringIO import StringIO
import os
import sys

model_prune_master_phil = """
  resolution_factor = 1/4.
    .type = float
    .help = Map grid spacing (multiplied times d_min).
  sidechains = True
    .type = bool
    .help = Remove poor sidechains
  residues = True
    .type = bool
    .help = Remove entire residues in poor density
  min_backbone_2fofc = 0.8
    .type = float
    .help = Minimum 2mFo-DFc sigma level at C-alpha to keep.  Residues with \
      C-alpha in density below this cutoff will be deleted.
  min_backbone_fofc = -3.0
    .type = float
    .help = Maximum mFo-DFc sigma level at C-alpha to keep.  Residues with \
      C-alpha in difference density below this cutoff will be deleted.
  min_sidechain_2fofc = 0.5
    .type = float
    .help = Minimum mean 2mFo-DFc sigma level for sidechain atoms to keep. \
      Residues with sidechains below this level will be truncated.
  max_sidechain_fofc = -2.8
    .type = float
    .help = Maximum mean 2mFo-DFc sigma level for sidechain atoms to keep. \
      Residues with sidechains below this level will be truncated.
  min_cc = 0.7
    .type = float
    .help = Minimum overall CC for entire residue to keep.
  min_cc_sidechain = 0.6
    .type = float
    .help = Minimum overall CC for sidechains to keep.
  min_fragment_size = 3
    .type = int
    .help = Minimum fragment size to keep.  Fragments smaller than this will \
      be deleted in the final step (based on the assumption that the adjacent \
      residues were already removed).  Set this to None to prevent fragment \
      filtering.
"""

master_phil = libtbx.phil.parse("""
include scope mmtbx.utils.cmdline_input_phil_str
prune {
  %s
}
output {
  file_name = None
    .type = path
}
""" % model_prune_master_phil, process_includes=True)

def id_str (chain, residue_group, atom_group) :
  return "%3s %s%4s%s" % (atom_group.resname, chain.id, residue_group.resseq,
    residue_group.icode)

class residue_summary (object) :
  def __init__ (self,
                chain_id,
                residue_group,
                atom_group,
                score,
                score_type="CC",
                map_type="2mFo-DFc",
                atoms_type="residue") :
    self.resname = atom_group.resname
    self.chain_id = chain_id
    self.resseq = residue_group.resseq
    self.icode = residue_group.icode
    self.score = score
    self.score_type = score_type
    self.atoms_type = atoms_type
    self.map_type = map_type

  def show (self, out=None) :
    if (out is None) : out = sys.stdout
    id_str = "%3s %s%4s%s" % (self.resname, self.chain_id, self.resseq,
      self.icode)
    if (self.score is not None) :
      print >> out, "%s : %s %s %s = %.2f" % (id_str, self.atoms_type,
        self.map_type, self.score_type, self.score)
    else :
      print >> out, "%s : not part of a continuous chain" % id_str

class prune_model (object) :
  def __init__ (self,
                f_map_coeffs,
                diff_map_coeffs,
                model_map_coeffs,
                pdb_hierarchy,
                params) :
    """
    Removes atoms with poor electron density, as judged by several sigma-level
    cutoffs and overall CC.  This is basically an attempt to apply the same
    visual intuition we use when editing in Coot (etc.).
    """
    # XXX or, viewed a different way, this is one giant hack.  need to make
    # the logic smarter!
    assert (len(pdb_hierarchy.models()) == 1)
    from mmtbx.real_space_correlation import set_details_level_and_radius
    from cctbx import maptbx
    from scitbx.array_family import flex
    self.params = params
    self.pdb_hierarchy = pdb_hierarchy
    self.unit_cell = f_map_coeffs.unit_cell()
    f_map_fft = f_map_coeffs.fft_map(resolution_factor=params.resolution_factor)
    self.f_map = f_map_fft.apply_sigma_scaling().real_map()
    diff_map_fft = diff_map_coeffs.fft_map(
      resolution_factor=params.resolution_factor)
    self.diff_map = diff_map_fft.apply_sigma_scaling().real_map()
    model_map_fft = model_map_coeffs.fft_map(
      resolution_factor=params.resolution_factor)
    self.model_map = model_map_fft.apply_sigma_scaling().real_map()
    atom_detail, residue_detail, self.atom_radius = \
      set_details_level_and_radius(
        details_level="automatic",
        d_min=f_map_coeffs.d_min(),
        atom_radius=None)

  def get_map_stats_for_atoms (self, atoms) :
    from cctbx import maptbx
    from scitbx.array_family import flex
    sites_cart = flex.vec3_double()
    sites_cart_nonH = flex.vec3_double()
    for atom in atoms :
      sites_cart.append(atom.xyz)
      if (not atom.element.strip() in ["H","D"]) :
        sites_cart_nonH.append(atom.xyz)
    if (len(sites_cart_nonH) == 0) :
      return None
    sel = maptbx.grid_indices_around_sites(
      unit_cell=self.unit_cell,
      fft_n_real=self.f_map.focus(),
      fft_m_real=self.f_map.all(),
      sites_cart=sites_cart,
      site_radii=get_atom_radii(atoms, self.atom_radius))
    f_map_sel = self.f_map.select(sel)
    model_map_sel = self.model_map.select(sel)
    diff_map_sel = self.diff_map.select(sel)
    mean_f_value = flex.mean(f_map_sel.as_1d())
    mean_diff_value = flex.mean(diff_map_sel.as_1d())
    cc = flex.linear_correlation(x=f_map_sel,
      y=model_map_sel).coefficient()
    return group_args(cc=cc,
      mean_2fofc=mean_f_value,
      mean_fofc=mean_diff_value)

  def process_residues (self, out=None) :
    if (out is None) :
      out = sys.stdout
    n_res_removed = 0
    n_sc_removed = 0
    n_res_protein = 0
    pruned = []
    make_header("Pruning residues and sidechains", out=out)
    for chain in self.pdb_hierarchy.models()[0].chains() :
      residue_id_hash = {}
      removed_resseqs = []
      if (len(chain.conformers()) > 1) :
        print >> out, "WARNING: chain '%s' has multiple conformers" % chain.id
      first_conf = chain.conformers()[0]
      if (not first_conf.is_protein()) :
        continue
      for j_seq, residue_group in enumerate(chain.residue_groups()) :
        n_res_protein += 1
        residue_id_hash[residue_group.resid()] = j_seq
        for atom_group in residue_group.atom_groups() :
          ag_id_str = id_str(chain, residue_group, atom_group)
          remove_atom_group = False
          sidechain_atoms = []
          backbone_atoms = []
          for atom in atom_group.atoms() :
            if (atom.name.strip() in ["N", "O", "C", "H", "CA", "CB"]) :
              backbone_atoms.append(atom)
            else :
              sidechain_atoms.append(atom)
          if (len(backbone_atoms) > 0) and (self.params.residues) :
            mc_stats = self.get_map_stats_for_atoms(backbone_atoms)
            if (mc_stats.mean_2fofc < self.params.min_backbone_2fofc) :
              pruned.append(residue_summary(
                chain_id=chain.id,
                residue_group=residue_group,
                atom_group=atom_group,
                score=mc_stats.mean_2fofc,
                score_type="sigma",
                atoms_type="C-alpha"))
              remove_atom_group = True
            elif (mc_stats.mean_fofc < self.params.min_backbone_fofc) :
              pruned.append(residue_summary(
                chain_id=chain.id,
                residue_group=residue_group,
                atom_group=atom_group,
                score=mc_stats.mean_fofc,
                score_type="sigma",
                atoms_type="C-alpha"))
              remove_atom_group = True
          # map values look okay - now check overall CC
          if (not remove_atom_group) :
            res_stats = self.get_map_stats_for_atoms(atom_group.atoms())
            if (res_stats.cc < self.params.min_cc) and (self.params.residues) :
              pruned.append(residue_summary(
                chain_id=chain.id,
                residue_group=residue_group,
                atom_group=atom_group,
                score=cc))
              remove_atom_group = True
            elif (len(sidechain_atoms) > 0) and (self.params.sidechains) :
              # overall CC is acceptable - now look at sidechain alone
              remove_sidechain = False
              sc_stats = self.get_map_stats_for_atoms(sidechain_atoms)
              if (sc_stats.cc < self.params.min_cc_sidechain) :
                pruned.append(residue_summary(
                  chain_id=chain.id,
                  residue_group=residue_group,
                  atom_group=atom_group,
                  score=sc_stats.cc,
                  atoms_type="sidechain"))
                remove_sidechain = True
              else :
                if (sc_stats.mean_2fofc < self.params.min_sidechain_2fofc) :
                  pruned.append(residue_summary(
                    chain_id=chain.id,
                    residue_group=residue_group,
                    atom_group=atom_group,
                    score=sc_stats.mean_2fofc,
                    score_type="sigma",
                    atoms_type="sidechain"))
                  remove_sidechain = True
                elif (sc_stats.mean_fofc < self.params.max_sidechain_fofc) :
                  pruned.append(residue_summary(
                    chain_id=chain.id,
                    residue_group=residue_group,
                    atom_group=atom_group,
                    score=sc_stats.mean_fofc,
                    score_type="sigma",
                    atoms_type="sidechain",
                    map_type="mFo-Dfc"))
                  remove_sidechain = True
              if (remove_sidechain) :
                assert (self.params.sidechains)
                for atom in sidechain_atoms :
                  atom_group.remove_atom(atom)
                n_sc_removed += 1
          if (remove_atom_group) :
            assert (self.params.residues)
            residue_group.remove_atom_group(atom_group)
        if (len(residue_group.atom_groups()) == 0) :
          chain.remove_residue_group(residue_group)
          n_res_removed += 1
          removed_resseqs.append(residue_group.resseq_as_int())
      # Final pass: remove lone single/pair residues
      if ((self.params.residues) and
          (self.params.min_fragment_size is not None)) :
        n_rg = len(chain.residue_groups())
        for j_seq, residue_group in enumerate(chain.residue_groups()) :
          if (residue_group.icode.strip() != "") :
            continue
          resseq = residue_group.resseq_as_int()
          remove = False
          if (resseq - 1 in removed_resseqs) or (j_seq == 0) :
            print "candidate:", resseq
            for k in range(1, self.params.min_fragment_size+1) :
              if (resseq + k in removed_resseqs) :
                remove = True
                break
              elif ((j_seq + k) >= len(chain.residue_groups())) :
                remove = True
                break
          if (remove) :
            pruned.append(residue_summary(
              chain_id=chain.id,
              residue_group=residue_group,
              atom_group=atom_group,
              score=None))
            chain.remove_residue_group(residue_group)
            removed_resseqs.append(resseq)
            n_res_removed += 1
    for outlier in pruned :
      outlier.show(out)
    print >> out, "Removed %d residues and %d sidechains" % (n_res_removed,
      n_sc_removed)
    return group_args(
      n_res_protein=n_res_protein,
      n_res_removed=n_res_removed,
      n_sc_removed=n_sc_removed,
      outliers=pruned)

def get_atom_radii (atoms, atom_radius) :
  from scitbx.array_family import flex
  radii = flex.double([atom_radius] * len(atoms))
  for i_seq, atom in enumerate(atoms) :
    if (atom.element.strip().upper() in ["H", "D"]) :
      radii[i_seq] = 1.0
  return radii

def run_post_refinement (
    pdb_file,
    map_coeffs_file,
    output_file=None,
    params=None,
    f_map_label="2FOFCWT,PH2FOFCWT",
    diff_map_label="FOFCWT,PHFOFCWT",
    model_map_label="F-model,PHIF-model",
    write_model=True,
    out=None) :
  if (out is None) : out = sys.stdout
  if (params is None) :
    params = master_phil.fetch().extract().prune
  from iotbx import file_reader
  pdb_in = file_reader.any_file(pdb_file, force_type="pdb")
  pdb_in.assert_file_type("pdb")
  pdb_hierarchy = pdb_in.file_object.construct_hierarchy()
  pdb_hierarchy.atoms().reset_i_seq()
  # XXX this probably shouldn't be necessary
  pdb_hierarchy.atoms().set_chemical_element_simple_if_necessary()
  mtz_in = file_reader.any_file(map_coeffs_file, force_type="hkl")
  mtz_in.assert_file_type("hkl")
  f_map_coeffs = diff_map_coeffs = model_map_coeffs = None
  for array in mtz_in.file_server.miller_arrays :
    labels = array.info().label_string()
    if (labels == f_map_label) :
      f_map_coeffs = array
    elif (labels == diff_map_label) :
      diff_map_coeffs = array
    elif (labels == model_map_label) :
      model_map_coeffs = array
  if (f_map_coeffs is None) :
    raise RuntimeError("2mFo-DFc map not found (expected labels %s)." %
      f_map_label)
  elif (diff_map_coeffs is None) :
    raise RuntimeError("mFo-DFc map not found (expected labels %s)." %
      diff_map_label)
  elif (model_map_coeffs is None) :
    raise RuntimeError("Fc map not found (expected labels %s)." %
      model_map_label)
  result = prune_model(
    f_map_coeffs=f_map_coeffs,
    diff_map_coeffs=diff_map_coeffs,
    model_map_coeffs=model_map_coeffs,
    pdb_hierarchy=pdb_hierarchy,
    params=params).process_residues(out=out)
  if (write_model) :
    if (output_file is None) :
      base_name = os.path.basename(pdb_file)
      output_file = os.path.splitext(base_name)[0] + "_pruned.pdb"
    f = open(output_file, "w")
    f.write("%s\n" % "\n".join(pdb_in.file_object.crystallographic_section()))
    f.write(pdb_hierarchy.as_pdb_string())
    f.close()
    result.output_file = output_file
  return result

def run (args, out=None) :
  if (out is None) : out = sys.stdout
  if (len(args) == 0) :
    phil_out = StringIO()
    master_phil.show(out=phil_out, attributes_level=1)
    raise Usage("""\
mmtbx.prune_model model.pdb data.mtz [options...]

Filters protein residues based on CC to 2mFo-DFc map and absolute
(sigma-scaled) values in 2mFo-DFc and mFo-DFc maps.  For fast automatic
correction of MR solutions after initial refinement (ideally with rotamer
correction) to remove spurious loops and sidechains.

Full parameters:
%s
""" % phil_out.getvalue())
  from mmtbx.utils import cmdline_load_pdb_and_data
  import iotbx.pdb
  cmdline = cmdline_load_pdb_and_data(
    args=args,
    master_phil=master_phil,
    out=out,
    process_pdb_file=False,
    create_fmodel=True)
  params = cmdline.params
  fmodel = cmdline.fmodel
  if (params.output.file_name is None) :
    base_name = os.path.basename(params.input.pdb.file_name[0])
    params.output.file_name = os.path.splitext(base_name)[0] + "_pruned.pdb"
  log_file = os.path.splitext(os.path.basename(params.output.file_name))[0] + \
    ".log"
  log = open(log_file, "w")
  out2 = multi_out()
  out2.register("out", out)
  out2.register("log", log)
  map_helper = fmodel.electron_density_map()
  f_map_coeffs = map_helper.map_coefficients(map_type="2mFo-DFc")
  diff_map_coeffs = map_helper.map_coefficients(map_type="mFo-DFc")
  model_map_coeffs = map_helper.map_coefficients(map_type="Fc")
  result = prune_model(
    f_map_coeffs=f_map_coeffs,
    diff_map_coeffs=diff_map_coeffs,
    model_map_coeffs=model_map_coeffs,
    pdb_hierarchy=cmdline.pdb_hierarchy,
    params=params.prune).process_residues(out=out2)
  f = open(params.output.file_name, "w")
  cryst1 = iotbx.pdb.format_cryst1_record(fmodel.xray_structure)
  f.write("REMARK edited by mmtbx.prune_model\n")
  f.write("%s\n" % cryst1)
  f.write(cmdline.pdb_hierarchy.as_pdb_string())
  f.close()
  log.close()
  print >> out, "Wrote %s" % params.output.file_name
  return params.output.file_name

if (__name__ == "__main__") :
  run(sys.argv[1:])
