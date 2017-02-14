from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.map_model_cc
# LIBTBX_SET_DISPATCHER_NAME phenix.model_map_cc

from scitbx.array_family import flex
import sys
import iotbx.pdb
from libtbx.utils import Sorry
import mmtbx.utils
import mmtbx.maps.correlation
from cctbx import maptbx
from cctbx.maptbx import resolution_from_map_and_model

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

def run(args, log=sys.stdout):
  print >> log, "-"*79
  print >> log, legend
  print >> log, "-"*79
  inputs = mmtbx.utils.process_command_line_args(args = args,
    master_params = master_params())
  params = inputs.params.extract()
  # model
  broadcast(m="Input PDB:", log=log)
  file_names = inputs.pdb_file_names
  if(len(file_names) != 1): raise Sorry("PDB file has to given.")
  pi = iotbx.pdb.input(file_name = file_names[0])
  h = pi.construct_hierarchy()
  xrs = h.extract_xray_structure(crystal_symmetry=inputs.crystal_symmetry)
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
  # Compute FSC(map, model)
  broadcast(m="Model-map FSC:", log=log)
  mmtbx.maps.correlation.fsc_model_map(
    xray_structure=xrs, map=map_data, d_min=d_min, log=log)
  #
  # various CC
  cc_calculator = mmtbx.maps.correlation.from_map_and_xray_structure_or_fmodel(
    xray_structure = xrs,
    map_data       = map_data,
    d_min          = d_min)
  broadcast(m="Map-model CC:", log=log)
  # per residue
  print >> log, "Per residue:"
  for rg in h.residue_groups():
    cc = cc_calculator.cc(selection=rg.atoms().extract_i_seq())
    print >> log, "  chain id: %s resid %s: %6.4f"%(
      rg.parent().id, rg.resid(), cc)
  # per chain
  print >> log, "Per chain:"
  for chain in h.chains():
    print >> log, "  chain %s: %6.4f"%(chain.id, cc_calculator.cc(
      selection=chain.atoms().extract_i_seq()))
  print >> log, "Overall:"
  # entire box
  print >> log, "         box: %6.4f"%cc_calculator.cc()
  # all atoms
  print >> log, "around atoms: %6.4f"%cc_calculator.cc(
    selection=flex.bool(xrs.scatterers().size(),True))
  #

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
