from __future__ import division
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
  scattering_table = wk1995  it1992  *n_gaussian  neutron electron
    .type = choice
    .help = Scattering table (X-ray, neutron or electron)
  atom_radius = None
    .type = float
    .help = Atom radius for masking. If undefined then calculated automatically
  compute_cc_per_chain = True
    .type = bool
    .help = Compute local model-map CC for each chain
  compute_cc_per_residue = True
    .type = bool
    .help = Compute local model-map CC for each residue
  compute_fsc = True
    .type = bool
    .help = Compute FSC
  compute_cc_mask = True
    .type = bool
  compute_cc_volume = True
    .type = bool
  compute_cc_peaks = True
    .type = bool
}
"""

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=False)

class map_model_cc(object):
  def __init__(self, map_data, pdb_hierarchy, crystal_symmetry, params):
    adopt_init_args(self, locals())
    self.five_cc_result = None
    self.fsc = None
    self.cc_per_chain = []
    self.cc_per_residue = []

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
    self.five_cc_result = mmtbx.maps.correlation.five_cc(
      map               = map_data,
      xray_structure    = xrs,
      d_min             = self.params.resolution,
      compute_cc_mask   = self.params.compute_cc_mask,
      compute_cc_volume = self.params.compute_cc_volume,
      compute_cc_peaks  = self.params.compute_cc_peaks)
    # Atom radius
    self.atom_radius = mtriage.get_atom_radius(
      xray_structure   = xrs,
      d_min            = self.params.resolution,
      map_data         = map_data,
      crystal_symmetry = self.crystal_symmetry,
      radius           = self.params.atom_radius)
    # Model-map FSC
    if(self.params.compute_fsc):
      mtriage_params = mtriage.master_params().extract()
      mtriage_params.scattering_table = self.params.scattering_table
      mtriage_params.compute_d_model=False
      mtriage_params.compute_d99=False
      mtriage_params.radius_smooth = self.atom_radius
      self.fsc = mtriage.mtriage(
        map_data         = self.map_data,
        pdb_hierarchy    = self.pdb_hierarchy,
        crystal_symmetry = self.crystal_symmetry,
        params           = mtriage_params).run().get_results().fsc_curve_model.fsc
    #
    cc_calculator = mmtbx.maps.correlation.from_map_and_xray_structure_or_fmodel(
      xray_structure = xrs,
      map_data       = map_data,
      d_min          = self.params.resolution)
    #
    def get_common_data(atoms, atom_radius):
      sel = atoms.extract_i_seq()
      return group_args(
        b_iso_mean = flex.mean(atoms.extract_b()),
        occ_mean   = flex.mean(atoms.extract_occ()),
        n_atoms    = atoms.size(),
        cc         = cc_calculator.cc(selection = sel, atom_radius = atom_radius),
        xyz_mean   = atoms.extract_xyz().mean())
    # CC per chain
    if(self.params.compute_cc_per_chain):
      for chain in self.pdb_hierarchy.chains():
        cd = get_common_data(atoms=chain.atoms(), atom_radius=self.atom_radius)
        self.cc_per_chain.append(group_args(
          chain_id   = chain.id,
          b_iso_mean = cd.b_iso_mean,
          occ_mean   = cd.occ_mean,
          n_atoms    = cd.n_atoms,
          cc         = cd.cc))
    # CC per residue
    if(self.params.compute_cc_per_residue):
      for rg in self.pdb_hierarchy.residue_groups():
        for conformer in rg.conformers():
          for residue in conformer.residues():
            cd = get_common_data(atoms=residue.atoms(), atom_radius=self.atom_radius)
            self.cc_per_residue.append(group_args(
              chain_id   = rg.parent().id,
              resname    = residue.resname,
              resseq     = residue.resseq,
              icode      = residue.icode,
              b_iso_mean = cd.b_iso_mean,
              occ_mean   = cd.occ_mean,
              n_atoms    = cd.n_atoms,
              cc         = cd.cc,
              xyz_mean   = cd.xyz_mean,
              residue    = residue)) # for compatibility with GUI

  def get_results(self):
    return group_args(
      cc_mask        = self.five_cc_result.cc_mask,
      cc_volume      = self.five_cc_result.cc_volume,
      cc_peaks       = self.five_cc_result.cc_peaks,
      cc_per_chain   = self.cc_per_chain,
      cc_per_residue = self.cc_per_residue,
      atom_radius    = self.atom_radius,
      fsc            = self.fsc)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
