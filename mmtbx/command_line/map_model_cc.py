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

legend = """phenix.map_model_cc or phenix.model_map_cc:
  Given PDB file and a map file calculate model-map coorelation.

How to run:
  phenix.map_model_cc model.pdb map.ccp4
  phenix.model_map_cc model.pdb map.ccp4

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
  xrs = pi.xray_structure_simple(crystal_symmetry=inputs.crystal_symmetry)
  xrs.scattering_type_registry(table = params.scattering_table)
  xrs.show_summary(f=log, prefix="  ")
  # map
  broadcast(m="Input map:", log=log)
  if(inputs.ccp4_map is None): raise Sorry("Map file has to given.")
  inputs.ccp4_map.show_summary(prefix="  ")
  map_data = inputs.ccp4_map.data.as_double()
  # estimate resolution
  d_min = params.resolution
  if(d_min is None):
    broadcast(m="Map resolution estimate:", log=log)
    d_min = maptbx.resolution_from_map_and_model(
      map_data=map_data, xray_structure=xrs)
  print >> log, "  d_min: %6.4f"%d_min
  # various CC
  cc_calculator = mmtbx.maps.correlation.from_map_and_xray_structure_or_fmodel(
    xray_structure = xrs,
    map_data       = map_data,
    d_min          = d_min)
  broadcast(m="Map-model CC:", log=log)
  # entire box
  print >> log, "      box: %6.4f"%cc_calculator.cc()
  # all atoms
  print >> log, "      all: %6.4f"%cc_calculator.cc(
    selection=flex.bool(xrs.scatterers().size(),True))
  # per chain
  for chain in h.chains():
    print >> log, "  chain %s: %6.4f"%(chain.id, cc_calculator.cc(
      selection=chain.atoms().extract_i_seq()))
  # per residue
  for rg in h.residue_groups():
    cc = cc_calculator.cc(selection=rg.atoms().extract_i_seq())
    print >> log, "  resseq %s: %6.4f"%(rg.resseq, cc)
  #

  ####
  #### SPECIAL IDIOTIC CASE:
  ####
#  cc1 = []
#  for i in xrange(160):
#    f = flex.double()
#    cc1.append(f)
#  cc2 = []
#  for i in xrange(182):
#    f = flex.double()
#    cc2.append(f)
#  #
#  for chain in h.chains():
#    n_at = len(chain.atoms())
#    if(n_at==1243):  # 160 residues
#      cntr = 0
#      for r in chain.residues():
#        sel = r.atoms().extract_i_seq()
#        cc1[cntr].append(cc_calculator.cc(selection=sel))
#        cntr+=1
#    if(n_at==1466):  # 160 residues
#      cntr = 0
#      for r in chain.residues():
#        sel = r.atoms().extract_i_seq()
#        cc2[cntr].append(cc_calculator.cc(selection=sel))
#        cntr+=1
#  #
#  for i, ccs in enumerate(cc1):
#    ccs = ccs.set_selected(ccs<0, 0)
#    v1,v2,v3=ccs.min_max_mean().as_tuple()
#    v4=ccs.sample_standard_deviation()
#    print "%3d: %6.4f %6.4f %6.4f %6.4f"%(i+1, v1,v2,v3,v4)
#  print
#  for i, ccs in enumerate(cc2):
#    ccs = ccs.set_selected(ccs<0, 0)
#    v1,v2,v3=ccs.min_max_mean().as_tuple()
#    v4=ccs.sample_standard_deviation()
#    print "%3d: %6.4f %6.4f %6.4f %6.4f"%(i+1, v1,v2,v3,v4)


if (__name__ == "__main__"):
  run(args=sys.argv[1:])
