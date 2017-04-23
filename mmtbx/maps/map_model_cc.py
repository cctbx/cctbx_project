from __future__ import division
from scitbx.array_family import flex
import sys
from libtbx import group_args
from libtbx.utils import Sorry
import mmtbx.maps.correlation
from cctbx import maptbx
import iotbx.phil

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
}
"""

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=False)

class map_model_cc(object):
  def __init__(self, map_data, pdb_hierarchy, crystal_symmetry, params):
    # XXX No default value for params because it cannot work with default
    # resolution=None
    self.map_data = map_data
    self.pdb_hierarchy = pdb_hierarchy
    self.crystal_symmetry = crystal_symmetry
    self.params = params

  def validate(self):
    assert not None in [self.map_data, self.pdb_hierarchy, self.crystal_symmetry, self.params]
    if(self.params.resolution is None):
      raise Sorry("Resolution is required.")

  def run(self):
    # assert len(locals().keys()) == 4 # intentional done in validate()
    xrs = self.pdb_hierarchy.extract_xray_structure(
      crystal_symmetry=self.crystal_symmetry)
    xrs.scattering_type_registry(table = self.params.scattering_table)
    soi = maptbx.shift_origin_if_needed(map_data=self.map_data, xray_structure=xrs)
    map_data = soi.map_data
    xrs = soi.xray_structure
    self.five_cc_result = mmtbx.maps.correlation.five_cc(
      map            = map_data,
      xray_structure = xrs,
      d_min          = self.params.resolution)
    atom_radius = self.params.atom_radius
    if(atom_radius is None):
      atom_radius = self.five_cc_result.atom_radius
    self.fsc = mmtbx.maps.correlation.fsc_model_vs_map(
      xray_structure = xrs,
      map            = map_data,
      atom_radius    = atom_radius,
      d_min          = self.params.resolution)
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
        cc         = cc_calculator.cc(selection = sel, atom_radius = atom_radius))
    # CC per chain
    self.cc_per_chain = []
    for chain in self.pdb_hierarchy.chains():
      cd = get_common_data(atoms=chain.atoms(), atom_radius=atom_radius)
      self.cc_per_chain.append(group_args(
        chain_id   = chain.id,
        b_iso_mean = cd.b_iso_mean,
        occ_mean   = cd.occ_mean,
        n_atoms    = cd.n_atoms,
        cc         = cd.cc))
    # CC per residue
    self.cc_per_residue = []
    for rg in self.pdb_hierarchy.residue_groups():
      for conformer in rg.conformers():
        for residue in conformer.residues():
          cd = get_common_data(atoms=residue.atoms(), atom_radius=atom_radius)
          self.cc_per_residue.append(group_args(
            chain_id   = chain.id,
            resname    = residue.resname,
            resseq     = residue.resseq,
            icode      = residue.icode,
            b_iso_mean = cd.b_iso_mean,
            occ_mean   = cd.occ_mean,
            n_atoms    = cd.n_atoms,
            cc         = cd.cc))
  def get_results(self):
    return group_args(
      resolution     = self.params.resolution, # Not clear why input parameter is returned as result
      cc_mask        = self.five_cc_result.cc_mask,
      cc_volume      = self.five_cc_result.cc_volume,
      cc_peaks       = self.five_cc_result.cc_peaks,
      cc_per_chain   = self.cc_per_chain,
      cc_per_residue = self.cc_per_residue,
      fsc            = self.fsc)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
