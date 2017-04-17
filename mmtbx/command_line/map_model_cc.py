from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.map_model_cc
# LIBTBX_SET_DISPATCHER_NAME phenix.model_map_cc

from scitbx.array_family import flex
import sys
import iotbx.pdb
from libtbx import group_args
from libtbx.utils import Sorry
import mmtbx.utils
import mmtbx.maps.correlation

legend = """phenix.map_model_cc or phenix.model_map_cc:
  Given PDB file and a map file calculate model-map coorelation.

How to run:
  phenix.map_model_cc model.pdb map.ccp4 resolution=3
  phenix.model_map_cc model.pdb map.ccp4 resolution=3

Feedback:
  PAfonine@lbl.gov
  phenixbb@phenix-online.org"""

master_params_str = """
  map_file_name = None
    .type = str
  model_file_name = None
    .type = str
  resolution = None
    .type = float
  scattering_table = wk1995  it1992  *n_gaussian  neutron electron
    .type = choice
"""

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=False)

def broadcast(m, log):
  print >> log, "-"*79
  print >> log, m
  print >> log, "*"*len(m)

def run(args, pdb_hierarchy=None, crystal_symmetry=None, log=sys.stdout):
  print >> log, "-"*79
  print >> log, legend
  print >> log, "-"*79
  inputs = mmtbx.utils.process_command_line_args(args = args,
    master_params = master_params())
  params = inputs.params.extract()
  # model
  h = pdb_hierarchy
  if (h is None):
    broadcast(m="Input PDB:", log=log)
    file_names = inputs.pdb_file_names
    if(len(file_names) != 1): raise Sorry("PDB file has to given.")
    pi = iotbx.pdb.input(file_name = file_names[0])
    h = pi.construct_hierarchy()
  if (crystal_symmetry is None):
    crystal_symmetry = inputs.crystal_symmetry
  xrs = h.extract_xray_structure(crystal_symmetry=crystal_symmetry)
  xrs.scattering_type_registry(table = params.scattering_table)
  xrs.show_summary(f=log, prefix="  ")
  # map
  broadcast(m="Input map:", log=log)
  if(inputs.ccp4_map is None): raise Sorry("Map file has to given.")
  inputs.ccp4_map.show_summary(prefix="  ")
  map_data = inputs.ccp4_map.map_data()
  # shift origin if needed
  shift_needed = not \
    (map_data.focus_size_1d() > 0 and map_data.nd() == 3 and
     map_data.is_0_based())
  if(shift_needed):
    N = map_data.all()
    O=map_data.origin()
    map_data = map_data.shift_origin()
    # apply same shift to the model
    a,b,c = xrs.crystal_symmetry().unit_cell().parameters()[:3]
    sites_cart = xrs.sites_cart()
    sx,sy,sz = a/N[0]*O[0], b/N[1]*O[1], c/N[2]*O[2]
    sites_cart_shifted = sites_cart-\
      flex.vec3_double(sites_cart.size(), [sx,sy,sz])
    xrs.set_sites_cart(sites_cart_shifted)
  # estimate resolution
  d_min = params.resolution
  broadcast(m="Map resolution:", log=log)
  if(d_min is None):
    raise Sorry("Resolution is required.")
  print >> log, "  d_min:", d_min
  # Compute CC
  broadcast(m="Map-model CC (overall):", log=log)
  five_cc_result = mmtbx.maps.correlation.five_cc(map = map_data,
    xray_structure = xrs, d_min = d_min)
  print >> log, "  CC_mask  : %6.4f"%five_cc_result.cc_mask
  print >> log, "  CC_volume: %6.4f"%five_cc_result.cc_volume
  print >> log, "  CC_peaks : %6.4f"%five_cc_result.cc_peaks
  # Compute FSC(map, model)
  broadcast(m="Model-map FSC:", log=log)
  fsc = mmtbx.maps.correlation.fsc_model_vs_map(
    xray_structure = xrs,
    map            = map_data,
    atom_radius    = five_cc_result.atom_radius,
    d_min          = d_min)
  fsc.show(prefix="  ")
  # Local CC
  cc_calculator = mmtbx.maps.correlation.from_map_and_xray_structure_or_fmodel(
    xray_structure = xrs,
    map_data       = map_data,
    d_min          = d_min)
  broadcast(m="Map-model CC (local):", log=log)
  # per residue
  print >> log, "Per residue:"
  residue_results = list()
  for rg in h.residue_groups():
    cc = cc_calculator.cc(selection=rg.atoms().extract_i_seq())
    chain_id = rg.parent().id
    print >> log, "  chain id: %s resid %s: %6.4f"%(
      chain_id, rg.resid(), cc)
    # GUI output, use first conformer
    residue = rg.conformers()[0].residues()[0]
    residue_results.append( group_args(
      residue = residue,
      chain_id = chain_id,
      id_str = "%2s %1s %3s %4s %1s"%(
        chain_id, rg.conformers()[0].altloc, residue.resname, residue.resseq,
        residue.icode),
      cc = cc,
      map_value_1 = None,
      map_value_2 = None,
      b = flex.mean(residue.atoms().extract_b()),
      occupancy = flex.mean(residue.atoms().extract_occ()),
      n_atoms = residue.atoms().size()
    ) )
  # per chain
  print >> log, "Per chain:"
  for chain in h.chains():
    print >> log, "  chain %s: %6.4f"%(chain.id, cc_calculator.cc(
      selection=chain.atoms().extract_i_seq()))
  #
  return five_cc_result, residue_results

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
