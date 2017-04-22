from __future__ import division
from scitbx.array_family import flex
import sys
from libtbx import group_args
from libtbx.utils import Sorry
import mmtbx.maps.correlation
from cctbx import maptbx
import iotbx.phil

master_params_str = """
  map_file_name = None
    .type = str
    .help = Map file name
  model_file_name = None
    .type = str
    .help = Model file name
  resolution = None
    .type = float
    .help = Data (map) resolution
  scattering_table = wk1995  it1992  *n_gaussian  neutron electron
    .type = choice
    .help = Scattering table (X-ray, neutron or electron)
  atom_radius = None
    .type = float
    .help = Atom radius for masking. If undefined then calculated automatically
"""

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=False)

def run(map_data, pdb_hierarchy, crystal_symmetry, params=master_params()):
  assert len(locals().keys()) == 4 # intentional
  xrs = pdb_hierarchy.extract_xray_structure(
    crystal_symmetry=crystal_symmetry)
  xrs.scattering_type_registry(table = params.scattering_table)
  soi = maptbx.shift_origin_if_needed(map_data=map_data, xray_structure=xrs)
  map_data = soi.map_data
  xrs = soi.xray_structure
  resolution = params.resolution
  if(resolution is None):
    raise Sorry("Resolution is required.")
  five_cc_result = mmtbx.maps.correlation.five_cc(
    map            = map_data,
    xray_structure = xrs,
    d_min          = resolution)
  atom_radius = params.atom_radius
  if(atom_radius is None):
    atom_radius = five_cc_result.atom_radius
  fsc = mmtbx.maps.correlation.fsc_model_vs_map(
    xray_structure = xrs,
    map            = map_data,
    atom_radius    = atom_radius,
    d_min          = resolution)
  cc_calculator = mmtbx.maps.correlation.from_map_and_xray_structure_or_fmodel(
    xray_structure = xrs,
    map_data       = map_data,
    d_min          = resolution)
  #
  def get_common_data(atoms, atom_radius):
    sel = atoms.extract_i_seq()
    return group_args(
      b_iso_mean = flex.mean(atoms.extract_b()),
      occ_mean   = flex.mean(atoms.extract_occ()),
      n_atoms    = atoms.size(),
      cc         = cc_calculator.cc(selection = sel, atom_radius = atom_radius))
  # CC per chain
  cc_per_chain = []
  for chain in pdb_hierarchy.chains():
    cd = get_common_data(atoms=chain.atoms(), atom_radius=atom_radius)
    cc_per_chain.append(group_args(
      chain_id   = chain.id,
      b_iso_mean = cd.b_iso_mean,
      occ_mean   = cd.occ_mean,
      n_atoms    = cd.n_atoms,
      cc         = cd.cc))
  # CC per residue
  cc_per_residue = []
  for rg in pdb_hierarchy.residue_groups():
    for conformer in rg.conformers():
      for residue in conformer.residues():
        cd = get_common_data(atoms=residue.atoms(), atom_radius=atom_radius)
        cc_per_chain.append(group_args(
          chain_id   = chain.id,
          resname    = residue.resname,
          resseq     = residue.resseq,
          icode      = residue.icode,
          b_iso_mean = cd.b_iso_mean,
          occ_mean   = cd.occ_mean,
          n_atoms    = cd.n_atoms,
          cc         = cd.cc))
  #
  return group_args(
    resolution     = resolution,
    cc_mask        = five_cc_result.cc_mask,
    cc_volume      = five_cc_result.cc_volume,
    cc_peaks       = five_cc_result.cc_peaks,
    cc_per_chain   = cc_per_chain,
    cc_per_residue = cc_per_residue,
    fsc            = fsc)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
