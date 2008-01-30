import sys, math
from cctbx.array_family import flex
import iotbx.pdb
import mmtbx.f_model
from cctbx import miller
from mmtbx import find_peaks
from mmtbx import monomer_library
import mmtbx.monomer_library.pdb_interpretation
import mmtbx.monomer_library.server
from mmtbx import utils

class h_peak(object):
  def __init__(self, site_frac,
                     height,
                     dist,
                     site_frac_o,
                     atom_attribute_o):
    self.site_frac        = site_frac
    self.height           = height
    self.dist             = dist
    self.site_frac_o      = site_frac_o
    self.atom_attribute_o = atom_attribute_o

def angle(o,h1,h2):
  result = None
  a = h1[0]-o[0], h1[1]-o[1], h1[2]-o[2]
  b = h2[0]-o[0], h2[1]-o[1], h2[2]-o[2]
  ab = a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
  a_abs = math.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
  b_abs = math.sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2])
  a_abs_b_abs = a_abs * b_abs
  if(a_abs_b_abs != 0):
    result = math.acos(ab/(a_abs_b_abs))*180./math.pi
  return result

def find_hydrogens(fmodel,
                   atom_attributes_list,
                   map_type,
                   map_cutoff,
                   dist_min,
                   dist_max,
                   peak_search_level):
  fp_params = find_peaks.master_params.extract()
  fp_params.peak_search.peak_search_level = peak_search_level
  fp_params.map_next_to_model.min_model_peak_dist = dist_min
  fp_params.map_next_to_model.min_peak_peak_dist = 1.0
  fp_params.peak_search.min_cross_distance = 1.0
  fp_params.map_next_to_model.use_hydrogens = True
  fp_params.map_next_to_model.max_model_peak_dist = dist_max
  fp_manager = find_peaks.manager(fmodel     = fmodel,
                                  map_type   = map_type,
                                  map_cutoff = map_cutoff,
                                  params     = fp_params,
                                  log        = None)
  fp_manager.show_mapped(atom_attributes_list = atom_attributes_list)
  return fp_manager.peaks_mapped()

def extract_hoh_peaks(peaks, atom_attributes_list, xray_structure):
  scatterers = xray_structure.scatterers()
  assert scatterers.size() == len(atom_attributes_list)
  assert peaks.sites.size() == peaks.heights.size()
  assert peaks.heights.size() == peaks.iseqs_of_closest_atoms.size()
  perm = flex.sort_permutation(peaks.iseqs_of_closest_atoms)
  sites = peaks.sites.select(perm)
  heights = peaks.heights.select(perm)
  iseqs_of_closest_atoms = peaks.iseqs_of_closest_atoms.select(perm)

  print "\n\nWater H:\n"
  for s, h, i_seq in zip(sites, heights, iseqs_of_closest_atoms):
    aa = atom_attributes_list[i_seq]
    if(aa.resName in ["HOH","SOL","SOLV","WAT","DOD","TIP3"]):
      d = xray_structure.unit_cell().distance(s, scatterers[i_seq].site)
      print "peak= %8.3f closest distance to %s = %8.3f"%(
        h, aa.pdb_format(), d)
  ####
  print
  top = {}
  unit_cell = xray_structure.unit_cell()
  for s, h, i_seq in zip(sites, heights, iseqs_of_closest_atoms):
    aa = atom_attributes_list[i_seq]
    if(aa.resName in ["HOH","SOL","SOLV","WAT","DOD","TIP3"]):
      hp = h_peak(
        site_frac        = s,
        height           = h,
        dist             = unit_cell.distance(s, scatterers[i_seq].site),
        site_frac_o      = scatterers[i_seq].site,
        atom_attribute_o = aa)
      top.setdefault(i_seq,[]).append(hp)
  for key in top.keys():
    print atom_attributes_list[int(key)].pdb_format()
    sz = len(top[key])
    for j, jhp in enumerate(top[key]):
      if(sz == 1):
        print "  ", jhp.height, jhp.dist
      else:
        angles = []
        for k, khp in enumerate(top[key]):
          if(j != k):
            assert jhp.site_frac_o == khp.site_frac_o
            angl = angle(o=unit_cell.orthogonalize(jhp.site_frac_o),
                                h1=unit_cell.orthogonalize(jhp.site_frac),
                                h2=unit_cell.orthogonalize(khp.site_frac))
            angles.append(angl)
            print "  ", jhp.height, jhp.dist, angl


#def add_hydrogens():

#def exercise():
#  xrs_exact = iotbx.pdb.input(file_name = "m.pdb").xray_structure_simple()
#  xrs_part = iotbx.pdb.input(file_name = "m_omitH.pdb").xray_structure_simple()
#  miller_set = miller.build_set(
#    crystal_symmetry = xrs_exact.crystal_symmetry(),
#    anomalous_flag   = False,
#    d_min            = 1.0)
#  f_obs = abs(miller_set.structure_factors_from_scatterers(
#    xray_structure = xrs_exact,
#    algorithm      = "direct",
#    cos_sin_table  = False).f_calc())
#  sf_par = mmtbx.f_model.sf_and_grads_accuracy_params.extract()
#  sf_par.algorithm = "direct"
#  sf_par.cos_sin_table = False
#  flags = f_obs.array(data=flex.bool(f_obs.data().size(),False))
#  fmodel = mmtbx.f_model.manager(
#    xray_structure               = xrs_part,
#    sf_and_grads_accuracy_params = sf_par,
#    r_free_flags                 = flags,
#    target_name                  = "ls_wunit_k1",
#    f_obs                        = f_obs)
#  #
#  mon_lib_srv = monomer_library.server.server()
#  ener_lib = monomer_library.server.ener_lib()
#  processed_pdb_file = monomer_library.pdb_interpretation.process(
#    mon_lib_srv              = mon_lib_srv,
#    ener_lib                 = ener_lib,
#    file_name                = "m_omitH.pdb")
#  aal = processed_pdb_file.all_chain_proxies.stage_1.atom_attributes_list
#  #
#  peaks = find_hydrogens(fmodel = fmodel,
#                         atom_attributes_list = aal,
#                         map_type = "mFo-DFc",
#                         map_cutoff = 3,
#                         dist_min = 0.7,
#                         dist_max = 1.3)

def exercise_real(pdb_file_name = "ald.pdb",
                  cif_file_name = "ald.cif",
                  mtz_file_name = "ald.mtz",
                  pkl_file_name = "fmodel_ald.pickle"):
  if 1:
    cif_objects = None
    if(cif_file_name is not None):
      cif_objects = [].append((cif_file_name,
        mmtbx.monomer_library.server.read_cif(file_name = cif_file_name)))
    processed_pdb_files_srv = utils.process_pdb_file_srv(
      cif_objects = cif_objects)
    processed_pdb_file, pdb_inp = \
      processed_pdb_files_srv.process_pdb_files(
        pdb_file_names = [pdb_file_name])
    aal = processed_pdb_file.all_chain_proxies.stage_1.atom_attributes_list
    xrs = processed_pdb_file.xray_structure(show_summary = False)
    #
  if 0:
    from iotbx import reflection_file_reader
    reflection_file = reflection_file_reader.any_reflection_file(
      file_name = mtz_file_name)
    reflection_file_as_miller_arrays = \
      reflection_file.as_miller_arrays(merge_equivalents=True)
    for miller_array in reflection_file_as_miller_arrays:
      label=miller_array.info().labels[0]
      if label in ['I-obs']:
        f_obs = miller_array
      if label in ['R-free-flags']:
        flags = miller_array
    f_obs = f_obs.f_sq_as_f()
    flags = flags.array(data = flags.data() == 1)

    sf_par = mmtbx.f_model.sf_and_grads_accuracy_params.extract()
    sf_par.algorithm = "fft"
    sf_par.cos_sin_table = True
    sf_par.grid_resolution_factor=1/4.
    fmodel = mmtbx.f_model.manager(
      xray_structure               = xrs,
      sf_and_grads_accuracy_params = sf_par,
      r_free_flags                 = flags,
      target_name                  = "ls_wunit_k1",
      f_obs                        = f_obs)
    print fmodel.r_work(), fmodel.r_free()
    fmodel.update_solvent_and_scale()
    #
    from libtbx import easy_pickle
    easy_pickle.dump(pkl_file_name, (fmodel))
    #
  ###
  from libtbx import easy_pickle
  fmodel = easy_pickle.load(pkl_file_name)
  print fmodel.r_work(), fmodel.r_free()
  ##
  peaks = find_hydrogens(fmodel = fmodel,
                         atom_attributes_list = aal,
                         map_type = "mFo-DFc",
                         map_cutoff = 2,
                         dist_min = 0.7,
                         dist_max = 1.15,
                         peak_search_level = 3)
  extract_hoh_peaks(peaks = peaks, atom_attributes_list = aal, xray_structure = xrs)


if (__name__ == "__main__"):
  exercise_real()



