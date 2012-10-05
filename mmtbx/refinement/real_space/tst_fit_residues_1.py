from __future__ import division
import mmtbx.monomer_library.pdb_interpretation
import iotbx.mtz
from cctbx.array_family import flex
import time
from mmtbx import monomer_library
import mmtbx.refinement.real_space.fit_residues
import iotbx.pdb
from libtbx import group_args
import mmtbx.command_line.real_space_refine


def exercise(d_min = 2.0, resolution_factor = 0.25):
  # answer
  pdb_inp = iotbx.pdb.input(file_name="answer.pdb")#(source_info=None, lines=pdb_answer)
  #pdb_inp.write_pdb_file(file_name = "answer.pdb")
  xrs_answer = pdb_inp.xray_structure_simple()
  f_calc = xrs_answer.structure_factors(d_min = d_min).f_calc()
  fft_map = f_calc.fft_map(resolution_factor=resolution_factor)
  fft_map.apply_sigma_scaling()
  target_map = fft_map.real_map_unpadded()
  mtz_dataset = f_calc.as_mtz_dataset(column_root_label = "FCmap")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "answer.mtz")
  # poor
  mon_lib_srv = monomer_library.server.server()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv              = mon_lib_srv,
    ener_lib                 = monomer_library.server.ener_lib(),
    #raw_records              = flex.std_string(pdb_poor_str.splitlines()),
    file_name                = "poor.pdb",
    strict_conflict_handling = True,
    force_symmetry           = True,
    log                      = None)
  pdb_hierarchy_poor = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  xrs_poor = processed_pdb_file.xray_structure()
  sites_cart_poor = xrs_poor.sites_cart()
  #pdb_hierarchy_poor.write_pdb_file(file_name = "poor.pdb")
  dist = xrs_answer.mean_distance(other = xrs_poor)
  print "start:", dist
  ####
  target_map_object = group_args(
    data             = target_map,
    miller_array     = f_calc,
    crystal_gridding = fft_map)
  grm = mmtbx.command_line.real_space_refine.get_geometry_restraints_manager(
    processed_pdb_file = processed_pdb_file, xray_structure = xrs_poor)
  sm = mmtbx.refinement.real_space.structure_monitor(
    pdb_hierarchy               = pdb_hierarchy_poor,
    xray_structure              = xrs_poor,
    target_map_object           = target_map_object,
    geometry_restraints_manager = grm.geometry)
  sm.show(prefix="START")
  print
  sm.show_residues()
  ####

  for i in [1]:#[1,2,3]:
    result = mmtbx.refinement.real_space.fit_residues.manager(
      structure_monitor = sm,
      mon_lib_srv       = mon_lib_srv)
    print
    sm.show(prefix="MC %s"%str(i))
    print
    sm.show_residues(map_cc_all=0.9)
  #
  sm.pdb_hierarchy.write_pdb_file(file_name = "refined.pdb")
  dist = xrs_answer.mean_distance(other = sm.xray_structure)
  print "final:", dist
  #assert dist < 0.0035, dist
  #print
  #sm.update(xray_structure = xrs_poor)
  #sm.show()
  #print
  #sm.show_residues()


if(__name__ == "__main__"):
  t0 = time.time()
  exercise()
  print "Time: %6.4f"%(time.time()-t0)
