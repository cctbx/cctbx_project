from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
import sys
from libtbx import group_args
from libtbx.utils import Sorry
import mmtbx.maps.correlation
from cctbx import maptbx
import iotbx.phil
from libtbx import adopt_init_args
from mmtbx.maps import mtriage


master_params_str = """
map_model_cc {
  resolution = None
    .type = float
    .help = Data (map) resolution
    .short_caption = Resolution (optional in GUI)
    .expert_level=1
  scattering_table = wk1995  it1992 n_gaussian  neutron *electron
    .type = choice
    .help = Scattering table (X-ray, neutron or electron)
    .expert_level=0
  atom_radius = None
    .type = float
    .help = Atom radius for masking. If undefined then calculated automatically
    .short_caption = Atom radius (optional)
    .expert_level=0
  keep_map_calc = False
    .type = bool
    .help = Keep model-calculated map
    .expert_level=3
  wrapping = None
    .type = bool
    .help = You can specify whether your maps wrap around outside boundaries
    .short_caption = Wrapping (optional)
    .expert_level=3
  ignore_symmetry_conflicts = False
    .type = bool
    .help = You can ignore the symmetry information (CRYST1) from \
            coordinate files. This may be necessary if your model has been\
            placed in a box with box_map for example.
    .expert_level=2

  compute {
    cc_per_chain = True
      .type = bool
      .help = Compute local model-map CC for each chain
      .short_caption = CC by chain
    cc_per_residue = True
      .type = bool
      .help = Compute local model-map CC for each residue
      .short_caption = CC by residue
    cc_per_residue_group = False
      .type = bool
      .help = Compute local model-map CC for each residue group
      .short_caption = CC by residue group
    fsc = True
      .type = bool
      .help = Compute FSC
      .short_caption = FSC
    cc_mask = True
      .type = bool
    cc_volume = True
      .type = bool
    cc_peaks = True
      .type = bool
    cc_box = True
      .type = bool
    cc_image = False
      .type = bool
  }
}
"""

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=False)

class map_model_cc(object):
  def __init__(self, map_data, pdb_hierarchy, crystal_symmetry, params):
    adopt_init_args(self, locals())
    self.five_cc_result = None
    self.cc_per_chain = []
    self.cc_per_residue = []
    self.cc_main_chain = None
    self.cc_side_chain = None
    self.cc_per_residue_group = flex.double()

  def validate(self):
    assert not None in [self.map_data, self.pdb_hierarchy,
      self.crystal_symmetry, self.params]
    if(self.params.resolution is None):
      raise Sorry("Resolution is required.")

  def run(self):
    xrs = self.pdb_hierarchy.extract_xray_structure(
      crystal_symmetry=self.crystal_symmetry)
    xrs.scattering_type_registry(table = self.params.scattering_table)
    soin = maptbx.shift_origin_if_needed(
      map_data         = self.map_data,
      sites_cart       = xrs.sites_cart(),
      crystal_symmetry = self.crystal_symmetry)
    map_data = soin.map_data
    xrs.set_sites_cart(soin.sites_cart)
    self.five_cc = mmtbx.maps.correlation.five_cc(
      map               = map_data,
      xray_structure    = xrs,
      keep_map_calc     = self.params.keep_map_calc,
      d_min             = self.params.resolution,
      compute_cc_mask   = self.params.compute.cc_mask,
      compute_cc_volume = self.params.compute.cc_volume,
      compute_cc_peaks  = self.params.compute.cc_peaks,
      compute_cc_box    = self.params.compute.cc_box,
      compute_cc_image  = self.params.compute.cc_image)
    # Atom radius
    self.atom_radius = mtriage.get_atom_radius(
      xray_structure = xrs,
      resolution     = self.params.resolution,
      radius         = self.params.atom_radius)
    #
    cc_calculator = mmtbx.maps.correlation.from_map_and_xray_structure_or_fmodel(
      xray_structure = xrs,
      map_data       = map_data,
      d_min          = self.params.resolution)
    #
    def get_common_data(atoms, atom_radius):
      sel = atoms.extract_i_seq()
      cc = cc_calculator.cc(selection = sel, atom_radius = atom_radius)
      return group_args(
        b_iso_mean = flex.mean(atoms.extract_b()),
        occ_mean   = flex.mean(atoms.extract_occ()),
        n_atoms    = atoms.size(),
        cc         = cc,
        xyz_mean   = atoms.extract_xyz().mean())
    # CC per chain
    if(self.params.compute.cc_per_chain):
      for chain in self.pdb_hierarchy.chains():
        cd = get_common_data(atoms=chain.atoms(), atom_radius=self.atom_radius)
        self.cc_per_chain.append(group_args(
          chain_id   = chain.id,
          b_iso_mean = cd.b_iso_mean,
          occ_mean   = cd.occ_mean,
          n_atoms    = cd.n_atoms,
          cc         = cd.cc))
    # CC per residue
    if(self.params.compute.cc_per_residue):
      for rg in self.pdb_hierarchy.residue_groups():
        for conformer in rg.conformers():
          for residue in conformer.residues():
            cd = get_common_data(
              atoms       = residue.atoms(),
              atom_radius = self.atom_radius)
            self.cc_per_residue.append(group_args(
              model_id   = rg.parent().parent().id,
              chain_id   = rg.parent().id,
              resname    = residue.resname,
              resseq     = residue.resseq,
              icode      = residue.icode,
              b_iso_mean = cd.b_iso_mean,
              occ_mean   = cd.occ_mean,
              n_atoms    = cd.n_atoms,
              cc         = cd.cc,
              xyz_mean   = cd.xyz_mean))
    # CC per residue group
    if(self.params.compute.cc_per_residue_group):
      for rg in self.pdb_hierarchy.residue_groups():
        rg_cc = cc_calculator.cc(
          selection   = rg.atoms().extract_i_seq(),
          atom_radius = self.atom_radius)
        self.cc_per_residue_group.append(rg_cc)
    # Side chain
    sel_mc_str = "protein and (name C or name N or name CA or name O or name CB)"
    asc = self.pdb_hierarchy.atom_selection_cache()
    sel_mc = asc.selection(sel_mc_str)
    sel_sc = ~sel_mc
    if(sel_mc.count(True)>0):
      self.cc_main_chain = get_common_data(
        atoms       = self.pdb_hierarchy.select(sel_mc).atoms(),
        atom_radius = self.atom_radius)
    if(sel_sc.count(True)>0):
      self.cc_side_chain = get_common_data(
        atoms       = self.pdb_hierarchy.select(sel_sc).atoms(),
        atom_radius = self.atom_radius)

  def get_results(self):
    return group_args(
      cc_mask        = self.five_cc.result.cc_mask,
      cc_volume      = self.five_cc.result.cc_volume,
      cc_peaks       = self.five_cc.result.cc_peaks,
      cc_box         = self.five_cc.result.cc_box,
      map_calc       = self.five_cc.result.map_calc,
      cc_per_chain   = self.cc_per_chain,
      cc_per_residue = self.cc_per_residue,
      cc_per_residue_group = self.cc_per_residue_group,
      cc_main_chain  = self.cc_main_chain,
      cc_side_chain  = self.cc_side_chain,
      atom_radius    = self.atom_radius)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
