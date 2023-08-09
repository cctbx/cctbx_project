
# TODO trim sidechains one atom at a time

from __future__ import absolute_import, division, print_function
from libtbx.str_utils import make_header
from libtbx.utils import multi_out
from libtbx import group_args
import os
import sys
from six.moves import range

model_prune_master_phil = """
  resolution_factor = 1/4.
    .type = float
    .help = Map grid spacing (multiplied times d_min).
  sidechains = True
    .type = bool
    .help = Remove poor sidechains
  mainchain = False
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
  min_sidechain_2fofc = 0.6
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
  check_cgamma = True
    .type = bool
    .help = Check for poor density at the C-gamma atom for long sidechains.  \
      Useful in cases where the terminal atoms may have been misfit into \
      nearby density.
"""

def get_master_phil():
  from mmtbx.command_line import generate_master_phil_with_inputs
  return generate_master_phil_with_inputs(
    enable_automatic_twin_detection=True,
    phil_string="""
prune {
  %s
}
output {
  file_name = None
    .type = path
}
""" % model_prune_master_phil)

def id_str(chain, residue_group, atom_group):
  return "%3s %s%4s%s" % (atom_group.resname, chain.id, residue_group.resseq,
    residue_group.icode)

class residue_summary(object):
  def __init__(self,
                chain_id,
                residue_group,
                atom_group,
                score,
                score_type="CC",
                map_type="2mFo-DFc",
                atoms_type="residue"):
    self.resname = atom_group.resname
    self.chain_id = chain_id
    self.resseq = residue_group.resseq
    self.icode = residue_group.icode
    self.score = score
    self.score_type = score_type
    self.atoms_type = atoms_type
    self.map_type = map_type

  def show(self, out=None):
    if (out is None) : out = sys.stdout
    id_str = "%3s %s%4s%s" % (self.resname, self.chain_id, self.resseq,
      self.icode)
    if (self.score is not None):
      print("%s : %s %s %s = %.2f" % (id_str, self.atoms_type,
        self.map_type, self.score_type, self.score), file=out)
    else :
      print("%s : not part of a continuous chain" % id_str, file=out)

class prune_model(object):
  def __init__(self,
                f_map_coeffs,
                diff_map_coeffs,
                model_map_coeffs,
                pdb_hierarchy,
                params):
    """
    Removes atoms with poor electron density, as judged by several sigma-level
    cutoffs and overall CC.  This is basically an attempt to apply the same
    visual intuition we use when editing in Coot (etc.).
    """
    # XXX or, viewed a different way, this is one giant hack.  need to make
    # the logic smarter!
    assert (len(pdb_hierarchy.models()) == 1)
    from mmtbx.real_space_correlation import set_detail_level_and_radius
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
    detail, self.atom_radius = \
      set_detail_level_and_radius(
        detail="automatic",
        d_min=f_map_coeffs.d_min(),
        atom_radius=None)

  def get_map_stats_for_atoms(self, atoms):
    from cctbx import maptbx
    from scitbx.array_family import flex
    sites_cart = flex.vec3_double()
    sites_cart_nonH = flex.vec3_double()
    values_2fofc = flex.double()
    values_fofc = flex.double()
    for atom in atoms :
      sites_cart.append(atom.xyz)
      if (not atom.element.strip() in ["H","D"]) : #XXX trap: neutrons?
        sites_cart_nonH.append(atom.xyz)
        site_frac = self.unit_cell.fractionalize(atom.xyz)
        values_2fofc.append(self.f_map.eight_point_interpolation(site_frac))
        values_fofc.append(self.diff_map.eight_point_interpolation(site_frac))
    if (len(sites_cart_nonH) == 0):
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
    cc = flex.linear_correlation(x=f_map_sel, y=model_map_sel).coefficient()
    return group_args(cc=cc,
      mean_2fofc=flex.mean(values_2fofc),
      mean_fofc=flex.mean(values_fofc))

  def get_density_at_atom(self, atom):
    site_frac = self.unit_cell.fractionalize(site_cart=atom.xyz)
    two_fofc_value = self.f_map.eight_point_interpolation(site_frac)
    fofc_value = self.diff_map.eight_point_interpolation(site_frac)
    return group_args(two_fofc=two_fofc_value, fofc=fofc_value)

  def process_residues(self, out=None):
    if (out is None):
      out = sys.stdout
    n_res_removed = 0
    n_sc_removed = 0
    n_res_protein = 0
    pruned = []
    make_header("Pruning residues and sidechains", out=out)
    for chain in self.pdb_hierarchy.models()[0].chains():
      if (not chain.is_protein()):
        continue
      residue_id_hash = {}
      removed_resseqs = []
      if (len(chain.conformers()) > 1):
        print("WARNING: chain '%s' has multiple conformers" % chain.id, file=out)
      for j_seq, residue_group in enumerate(chain.residue_groups()):
        n_res_protein += 1
        residue_id_hash[residue_group.resid()] = j_seq
        for atom_group in residue_group.atom_groups():
          ag_id_str = id_str(chain, residue_group, atom_group)
          resname = atom_group.resname
          remove_atom_group = False
          sidechain_atoms = []
          backbone_atoms = []
          for atom in atom_group.atoms():
            if (atom.name.strip() in ["N", "O", "C", "H", "CA", "CB"]):
              backbone_atoms.append(atom)
            elif (not atom_group.resname in ["ALA", "GLY"]):
              sidechain_atoms.append(atom)
          if (len(backbone_atoms) > 0) and (self.params.mainchain):
            mc_stats = self.get_map_stats_for_atoms(backbone_atoms)
            if (mc_stats.mean_2fofc < self.params.min_backbone_2fofc):
              pruned.append(residue_summary(
                chain_id=chain.id,
                residue_group=residue_group,
                atom_group=atom_group,
                score=mc_stats.mean_2fofc,
                score_type="sigma",
                atoms_type="C-alpha"))
              remove_atom_group = True
            elif (mc_stats.mean_fofc < self.params.min_backbone_fofc):
              pruned.append(residue_summary(
                chain_id=chain.id,
                residue_group=residue_group,
                atom_group=atom_group,
                score=mc_stats.mean_fofc,
                score_type="sigma",
                atoms_type="C-alpha"))
              remove_atom_group = True
          # map values look okay - now check overall CC
          if (not remove_atom_group):
            res_stats = self.get_map_stats_for_atoms(atom_group.atoms())
            if (res_stats.cc < self.params.min_cc) and (self.params.mainchain):
              pruned.append(residue_summary(
                chain_id=chain.id,
                residue_group=residue_group,
                atom_group=atom_group,
                score=res_stats.cc))
              remove_atom_group = True
            elif (len(sidechain_atoms) > 0) and (self.params.sidechains):
              # overall CC is acceptable - now look at sidechain alone
              remove_sidechain = False
              sc_stats = self.get_map_stats_for_atoms(sidechain_atoms)
              if (sc_stats is None):
                continue
              if (sc_stats.cc < self.params.min_cc_sidechain):
                pruned.append(residue_summary(
                  chain_id=chain.id,
                  residue_group=residue_group,
                  atom_group=atom_group,
                  score=sc_stats.cc,
                  atoms_type="sidechain"))
                remove_sidechain = True
              else :
                if (sc_stats.mean_2fofc < self.params.min_sidechain_2fofc):
                  pruned.append(residue_summary(
                    chain_id=chain.id,
                    residue_group=residue_group,
                    atom_group=atom_group,
                    score=sc_stats.mean_2fofc,
                    score_type="sigma",
                    atoms_type="sidechain"))
                  remove_sidechain = True
                elif (sc_stats.mean_fofc < self.params.max_sidechain_fofc):
                  pruned.append(residue_summary(
                    chain_id=chain.id,
                    residue_group=residue_group,
                    atom_group=atom_group,
                    score=sc_stats.mean_fofc,
                    score_type="sigma",
                    atoms_type="sidechain",
                    map_type="mFo-Dfc"))
                  remove_sidechain = True
                if ((self.params.check_cgamma) and
                    (resname in ["ARG","LYS","TYR","TRP","PHE"])):
                  c_gamma = c_delta = None
                  for atom in atom_group.atoms():
                    if (atom.name.strip() == "CG"):
                      c_gamma = atom
                    elif (atom.name.strip() == "CD"):
                      c_delta = atom
                  if (c_gamma is not None):
                    map_values = self.get_density_at_atom(c_gamma)
                    # FIXME this is horribly subjective, but so is the logic
                    # I use for manual pruning...
                    if ((map_values.two_fofc < 0.8) or
                        ((map_values.two_fofc < 1.0) and
                         (map_values.fofc < -3.0))):
                      pruned.append(residue_summary(
                        chain_id=chain.id,
                        residue_group=residue_group,
                        atom_group=atom_group,
                        score=map_values.two_fofc,
                        score_type="sigma",
                        atoms_type="sidechain",
                        map_type="2mFo-Dfc"))
                      remove_sidechain = True
              if (remove_sidechain):
                assert (self.params.sidechains)
                for atom in sidechain_atoms :
                  atom_group.remove_atom(atom)
                n_sc_removed += 1
          if (remove_atom_group):
            assert (self.params.mainchain)
            residue_group.remove_atom_group(atom_group)
        if (len(residue_group.atom_groups()) == 0):
          chain.remove_residue_group(residue_group)
          n_res_removed += 1
          removed_resseqs.append(residue_group.resseq_as_int())
      # Final pass: remove lone single/pair residues
      if ((self.params.mainchain) and
          (self.params.min_fragment_size is not None)):
        n_rg = len(chain.residue_groups())
        for j_seq, residue_group in enumerate(chain.residue_groups()):
          if (residue_group.icode.strip() != ""):
            continue
          resseq = residue_group.resseq_as_int()
          remove = False
          if (resseq - 1 in removed_resseqs) or (j_seq == 0):
            print("candidate:", resseq)
            for k in range(1, self.params.min_fragment_size+1):
              if (resseq + k in removed_resseqs):
                remove = True
                break
              elif ((j_seq + k) >= len(chain.residue_groups())):
                remove = True
                break
          if (remove):
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
    print("Removed %d residues and %d sidechains" % (n_res_removed,
      n_sc_removed), file=out)
    return group_args(
      n_res_protein=n_res_protein,
      n_res_removed=n_res_removed,
      n_sc_removed=n_sc_removed,
      outliers=pruned)

def get_atom_radii(atoms, atom_radius):
  from scitbx.array_family import flex
  radii = flex.double([atom_radius] * len(atoms))
  for i_seq, atom in enumerate(atoms):
    if (atom.element.strip().upper() in ["H", "D"]):
      radii[i_seq] = 1.0
  return radii

def run_post_refinement(
    pdb_file,
    map_coeffs_file,
    output_file=None,
    params=None,
    f_map_label="2FOFCWT",
    diff_map_label="FOFCWT",
    model_map_label="F-model",
    write_model=True,
    out=None):
  if (out is None) : out = sys.stdout
  if (params is None):
    params = get_master_phil().fetch().extract().prune
  from iotbx import file_reader
  import iotbx.pdb
  pdb_in = iotbx.pdb.input(pdb_file)
  pdb_hierarchy = pdb_in.construct_hierarchy()
  pdb_hierarchy.atoms().reset_i_seq()
  # XXX this probably shouldn't be necessary
  pdb_hierarchy.atoms().set_chemical_element_simple_if_necessary()
  mtz_in = file_reader.any_file(map_coeffs_file, force_type="hkl")
  mtz_in.assert_file_type("hkl")
  f_map_coeffs = diff_map_coeffs = model_map_coeffs = None
  for array in mtz_in.file_server.miller_arrays :
    labels = array.info().labels
    if (labels[0] == f_map_label):
      f_map_coeffs = array
    elif (labels[0] == diff_map_label):
      diff_map_coeffs = array
    elif (labels[0] in [model_map_label, model_map_label + "(+)"]):
      model_map_coeffs = array.average_bijvoet_mates()
  if (f_map_coeffs is None):
    raise RuntimeError("2mFo-DFc map not found (expected labels %s)." %
      f_map_label)
  elif (diff_map_coeffs is None):
    raise RuntimeError("mFo-DFc map not found (expected labels %s)." %
      diff_map_label)
  elif (model_map_coeffs is None):
    raise RuntimeError("Fc map not found (expected labels %s)." %
      model_map_label)
  result = prune_model(
    f_map_coeffs=f_map_coeffs,
    diff_map_coeffs=diff_map_coeffs,
    model_map_coeffs=model_map_coeffs,
    pdb_hierarchy=pdb_hierarchy,
    params=params).process_residues(out=out)
  if (write_model):
    if (output_file is None):
      base_name = os.path.basename(pdb_file)
      output_file = os.path.splitext(base_name)[0] + "_pruned.pdb"
    f = open(output_file, "w")
    f.write("%s\n" % "\n".join(
      pdb_in.crystallographic_section()))
    f.write(pdb_hierarchy.as_pdb_string())
    f.close()
    result.output_file = output_file
  return result

def run(args, out=None):
  if (out is None) : out = sys.stdout
  usage_string = """\
mmtbx.prune_model model.pdb data.mtz [options...]

Filters protein residues based on CC to 2mFo-DFc map and absolute
(sigma-scaled) values in 2mFo-DFc and mFo-DFc maps.  For fast automatic
correction of MR solutions after initial refinement (ideally with rotamer
correction) to remove spurious loops and sidechains.
"""
  from mmtbx.command_line import load_model_and_data
  cmdline = load_model_and_data(
    args=args,
    master_phil=get_master_phil(),
    out=out,
    process_pdb_file=False,
    create_fmodel=True)
  params = cmdline.params
  fmodel = cmdline.fmodel
  if (params.output.file_name is None):
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
  f.write("REMARK edited by mmtbx.prune_model\n")
  f.write(cmdline.pdb_hierarchy.as_pdb_string(
    crystal_symmetry=fmodel.xray_structure))
  f.close()
  log.close()
  print("Wrote %s" % params.output.file_name, file=out)
  return params.output.file_name

if (__name__ == "__main__"):
  run(sys.argv[1:])
