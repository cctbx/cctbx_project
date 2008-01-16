# LIBTBX_SET_DISPATCHER_NAME phenix.cif_as_mtz

import sys, os, time
from cctbx.array_family import flex
from cctbx import miller
from libtbx import easy_run
from cctbx import crystal
from iotbx.option_parser import iotbx_option_parser
from libtbx.utils import Sorry
from iotbx import reflection_file_utils
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from cctbx.eltbx.neutron import neutron_news_1992_table
from cctbx import eltbx
import mmtbx.f_model
import mmtbx.bulk_solvent.bulk_solvent_and_scaling as bss
from iotbx.pdb import crystal_symmetry_from_pdb
from libtbx import smart_open
from libtbx.str_utils import show_string

"""
Notes on CIF (source: http://www.ccp4.ac.uk/html/mtz2various.html)
All reflections in the MTZ input file will be output to the CIF file. However,
there are ways to flag certain reflections with the data type _refln.status.
Observed reflections will be flagged with 'o'. Unobserved reflections, i.e.
those flagged as missing in the relevant amplitude or intensity column, will be
flagged as 'x'; these reflections will not be added to _reflns.number_obs. The
'free' reflections will be flagged as 'f'. The keyword FREEVAL can be used to
indicate this set. Systematically absent reflections are flagged with '-'.

If the RESO keyword is specified then reflections at higher or lower resolution
than the limits given, will be written with _refln.status 'h' or 'l'
respectively. The limits will be written to the CIF as the values of
_refine.ls_d_res_high and _refine.ls_d_res_low.

If EXCLUDE SIG is given then reflections for which F < <value>*sigma(F), and
which satisfy the resolution limits (if given), will be written with
_refln.status '<'. The value of _reflns.number_obs excludes all reflections
which do not satisfy the condition on sigma(F). All other sub-keywords of
EXCLUDE are ignored for CIF output.
NB: The translation of the RESOLUTION and EXCLUDE SIGP conditions to
_refln.status values does not imply that the the use of these conditions is
good crystallographic practice. Be prepared to justify why you have excluded
any data from your final refinement!
"""

keys_indices = ["index_h", "index_k", "index_l"]

keys_f_obs = ["F_meas_au",
              "F_meas",
              "F_meas_uncorrected"]
keys_sf_obs = ["F_meas_sigma_au",
               "F_meas_sigma",
               "F_sigma",
               "F_meas_siygma_au",
               "F_meas_au_sigma",
               "F_meas_sigma_uncorrected"]

keys_i_obs = ["intensity_meas",
              "F_squared_meas",
              "I_meas_au",
              "intensity_meas_au"]
keys_si_obs = ["intensity_sigma",
               "intensity_meas_sigma",
              "intensity__sigma",
              "I_meas_sigma_au",
              "F_squared_sigma",
              "F_squared_meas_sigma",
              "intensity_sigma_au",
              "intensity_sigm"]

keys_status = ["status", "R_free_flag", "statu", "observed_status", "F_status",
               "status_au"]

keys_other = ["phase_meas",
              "F.phase_calc",
              "pdbx_F_plus",
              "pdbx_F_plus_sigma",
              "pdbx_F_minus",
              "pdbx_F_minus_sigma",
              "pdbx_anom_difference",
              "pdbx_anom_difference_sigma",
              "pdbx_I_plus",
              "pdbx_I_plus_sigma",
              "pdbx_I_minus",
              "pdbx_I_minus_sigma",
              "pdbx_HL_A_iso",
              "pdbx_HL_B_iso",
              "pdbx_HL_C_iso",
              "pdbx_HL_D_iso",
              "pdbx_SAD_HL_A_iso",
              "pdbx_SAD_HL_B_iso",
              "pdbx_SAD_HL_C_iso",
              "pdbx_SAD_HL_D_iso",
              "pdbx_F_meas_plus_au",
              "pdbx_F_meas_plus_sigma_au",
              "pdbx_F_meas_minus_au",
              "pdbx_F_meas_minus_sigma_au",
              "Fc_au",
              "intensity_calc",
              "A_cal",
              "B_cal",
              "F-_meas_au",
              "F-_meas_sigma_au",
              "ccp4_I_plus",
              "ccp4_I_minus",
              "ccp4_I_plus_sigma",
              "ccp4_I_minus_sigma",
              "fiber_F_meas_au",
              "fiber_coordinate",
              "fiber_layer",
              "sgx_fmap",
              "pahse_calc",
              "A_cal_au",
              "B_cal_au",
              "ndb_anomalous_diff",
              "ndb_anomalous_diff_sigma",
              "pdbx_F_backtransform_au",
              "pdbx_phase_backtransform",
              "pdbx_anomalous_meas_au",
              "pdbx_anomalous_meas_sigma_au",
              "pdbx_SAD_phase_anom",
              "pdbx_SAD_phase_anom_sigma",
              "pdbx_phase_anom",
              "pdbx_phase_anom_sigma",
              "wavelength_di",
              "pdbx_phase_HLA",
              "pdbx_phase_HLB",
              "pdbx_phase_HLC",
              "pdbx_phase_HLD",
              "pdbx_fom_weighted_fmap",
              "ccp4_phase_anom",
              "ccp4_phase_anom_sigma",
              "ccp4_F_meas_minus_au",
              "ccp4_F_meas_minus_sigma_au",
              "ccp4_F_meas_plus_au",
              "ccp4_F_meas_plus_sigma_au",
              "ccp4_HL_A_iso",
              "ccp4_HL_B_iso",
              "ccp4_HL_C_iso",
              "ccp4_HL_D_iso",
              "pdbx_HLA",
              "pdbx_HLB",
              "pdbx_HLC",
              "pdbx_HLD",
              "phase_au",
              "DF_anomalous",
              "DF_anomalous_sigma",
              "d_spacing",
              "ccp4_SAD_phase_anom",
              "ccp4_SAD_phase_anom_sigma",
              "phase_meas_sigma",
              "F+_meas_au",
              "F+_meas_sigma_au",
              "F-_meas_au,"
              "F-_meas_sigma_au",
              "wavelength_id 1",
              "crystal_id 1",
              "pdbx_gsas_i100_meas",
              "ccp4_anomalous_diff",
              "ccp4_anomalous_diff_sigma",
              "phase_calc",
              "I_ano_meas",
              "I_ano_sigma",
              "wavelength_id",
              "crystal_id",
              "F_part_au",
              "weight",
              "F_calc",
              "phase_part",
              "scale_group_code",
              "fom",
              "ccp4_SAD_F_meas_plus_au",
              "ccp4_SAD_F_meas_plus_sigma_au",
              "ccp4_SAD_F_meas_minus_au",
              "ccp4_SAD_F_meas_minus_sigma_au",
              "ccp4_SAD_HL_A_iso",
              "ccp4_SAD_HL_B_iso",
              "ccp4_SAD_HL_C_iso",
              "ccp4_SAD_HL_D_iso",
              "F_calc_au",
              "A_calc",
              "B_calc",
              "A_calc_au",
              "B_calc_au",
              "resolution",
              "sint_over_lambda",
              "ebi_F_xplor_bulk_solvent_Calc",
              "ebi_Phase_xplor_bulk_solvent_Calc",
              "ebi_F_xplor_bulk_solvent_calc",
              "ebi_phase_xplor_bulk_solvent_calc",
              "Phase_Calc",
              "dfano",
              "pdbx_phase_cycle",
              "pdbx_sin_phase_calc",
              "pdbx_cos_phase_calc",
              "F_squared_calc"]
possible_keys = keys_indices + keys_f_obs + keys_sf_obs + keys_i_obs + \
                keys_si_obs + keys_status + keys_other

ccp4_range = range(-99,-1)+ range(2,250)
tw_flags_num = ["1","1.","1.0","1.00","1.000","1.0000","-1","-1.","-1.0","-1.00",
            "-1.000","0","0.0","0.","0.00","0.000","0.0000"]+ \
    [str(i) for i in ccp4_range] + ccp4_range + \
    [str("%.1f"%i) for i in ccp4_range]+[str("%.2f"%i) for i in ccp4_range]+\
    [str("%.3f"%i) for i in ccp4_range]+[str("%.1f"%i)[:-1] for i in ccp4_range]

cif_flags = ["x","o","f","-",">","<","h","l","O","F"]

tw_flags_all = tw_flags_num + cif_flags


def extract_keys(file_name, file_lines):
  lock_0 = False
  lock_1 = False
  keys = []
  for i_seq, st in enumerate(file_lines):
    st = st.strip()
    if(st.startswith("loop_")):
      if(i_seq < len(file_lines)):
        if(file_lines[i_seq+1].strip().startswith("_refln.") or
           file_lines[i_seq+1].strip().startswith("#_refln.")):
          lock_0 = True
    if(lock_0 and (st.startswith("_refln.") or st.startswith("#_refln."))):
      lock_1 = True
    if(lock_0 and lock_1 and not (st.startswith("_refln.") or
       st.startswith("#_refln."))):
      lock_0 = False
      lock_1 = False
    if(lock_0 and lock_1 and st[0] != "#"):
      keys.append(st)
    if(len(keys) > 0 and st.startswith("loop_")):
      print "Multiple data sets found => ignored: ", file_name
      return []
  if(len(keys) < 4):
    print "No keys found: ", file_name
    return []
  result = []
  if(len(keys) > 0):
    for i_seq, key in enumerate(keys):
      key = key.strip()
      key = key.replace("_refln.","")
      result.append(key)
  return result

class count_keys(object):
  def __init__(self, keys, file_name):
    self.keys = keys
    self.i_h = None
    self.i_k = None
    self.i_l = None
    self.i_fobs = None
    self.i_iobs = None
    self.i_sfobs = None
    self.i_siobs = None
    self.i_flag = None
    self.i_wavelength_id = None
    self.i_crystal_id = None
    n_h = 0
    n_k = 0
    n_l = 0
    n_fobs = 0
    n_iobs = 0
    n_flags = 0
    n_wavelength_id = 0
    n_crystal_id = 0
    for i_seq, key in enumerate(keys):
      if(key == "index_h"):
        self.i_h = i_seq
        n_h += 1
      elif(key == "index_k"):
        self.i_k = i_seq
        n_k += 1
      elif(key == "index_l"):
        self.i_l = i_seq
        n_l += 1
      elif(key in keys_f_obs):
        self.i_fobs = i_seq
        n_fobs += 1
      elif(key in keys_i_obs):
        self.i_iobs = i_seq
        n_iobs += 1
      elif(key in keys_sf_obs):
        self.i_sfobs = i_seq
      elif(key in keys_si_obs):
        self.i_siobs = i_seq
      elif(key in keys_status):
        self.i_flag = i_seq
        n_flags += 1
      elif(key in "wavelength_id"):
        self.i_wavelength_id = i_seq
        n_wavelength_id += 1
      elif(key in "crystal_id"):
        self.i_crystal_id = i_seq
        n_crystal_id += 1
      else:
        if(key not in possible_keys):
          proceed = False
          print "Unknown key= ", key, keys, file_name
          self.reset()
          break
    if(n_h*n_k*n_l != 1):
      print "No of multiple sets of Miller indices: ", keys, file_name
      self.reset()
    if(n_flags > 1):
      print "Multiple test sets: ", keys, file_name
      self.reset()
    if([self.i_fobs,self.i_iobs].count(None) == 2 or n_fobs > 1 or n_iobs > 1):
       print "Multiple or no data sets: ", keys, file_name
       self.reset()

  def reset(self):
    self.i_h, self.i_k, self.i_l = None, None, None

def get_array_of_r_free_flags(flags, crystal_symmetry, indices, file_name):
  flag_values = []
  guess_cif = 0
  cif_selection = None
  result = flex.int()
  for flag in flags:
    if(flag not in tw_flags_all):
      print "Unknown flag:", str(flag), file_name, " ... proceed anyway."
      return flex.int(), None
    if(flag not in flag_values):
      flag_values.append(flag)
  if(len(flag_values) == 0):
    print "No flag values:", file_name, " ... proceed anyway."
    return flex.int(), None
  if(len(flag_values) == 1): return flex.int(), None
  for item in flag_values:
    if(item in cif_flags): guess_cif += 1
  if(guess_cif > 0):
    cif_selection = flex.bool()
  for flag in flags:
    flag = flag.strip()
    if(guess_cif == 0):
      if(str(flag) in ["o","O"]): flag = 0
      result.append(int(float(flag)))
    if(guess_cif > 0):
      if(flag in ["f","F"]):
        result.append(1)
        cif_selection.append(True)
      elif(flag in ["o","O"]):
        result.append(0)
        cif_selection.append(True)
      elif(flag in ["-","x"]):
        result.append(0)
        cif_selection.append(False)
      elif(flag in tw_flags_num):
        result.append(int(float(flag)))
        cif_selection.append(True)
      else:
        result.append(0)
        cif_selection.append(True)
  flags_mi = None
  try:
    flags_mi = miller.set(
      crystal_symmetry = crystal_symmetry,
      indices          = indices).array(data = result)
  except Exception, e:
    print "Cannot setup miller array with flags: ", file_name, str(e), \
      " ... proceed anyway."
    return flex.int(), None
  if(flags_mi is not None):
    try:
      test_flag_value = reflection_file_utils.get_r_free_flags_scores(
        miller_arrays   = [flags_mi],
        test_flag_value = None).test_flag_values[0]
      if(test_flag_value is not None):
        result = (result == test_flag_value)
      else:
        print "Cannot score (reflection_file_utils.get_r_free_flags_scores):", \
          file_name, " ... proceed anyway."
        return flex.int(), None
    except Exception, e:
      print "Cannot determine flag value:", file_name, " ... proceed anyway.",\
        str(e)
      return flex.int(), None
  return result, cif_selection

class extract_data(object):
  def __init__(self, key_counter,
                     file_name,
                     file_lines,
                     wavelength_id,
                     crystal_id,
                     crystal_symmetry):
    self.indices = flex.miller_index()
    self.data = flex.double()
    self.sigmas = flex.double()
    self.flags = flex.std_string()
    self.file_name = file_name
    start_counter = 0
    start_flag = False
    for i_line, line in enumerate(file_lines):
      line_orig = line
      h_,k_,l_,data_,sigma_,flag_ = [None]*6
      line = line.strip()
      line = line.split()
      if len(line) == 0: continue
      if(len(key_counter.keys) != start_counter or not start_flag):
        if(len(line) == 1):
          if(line[0] == "loop_"): start_flag = True
          if(start_flag and (line[0].replace("_refln.","") in key_counter.keys)
             and line[0] != "loop_"):
            start_counter += 1
      if(len(key_counter.keys) == start_counter and start_flag):
        result_hkl = list(self.access_hkl(line=line, key_counter=key_counter))
        if(result_hkl.count(None) == 0):
          wavelength_id_and_crystal_id_ok = True
          if(wavelength_id is not None):
            if(int(line[key_counter.i_wavelength_id]) != wavelength_id):
              wavelength_id_and_crystal_id_ok = False
          if(crystal_id is not None):
            if(int(line[key_counter.i_crystal_id]) != crystal_id):
              wavelength_id_and_crystal_id_ok = False
          if(wavelength_id_and_crystal_id_ok):
            if(len(line) != len(key_counter.keys)):
              if(self.is_break(file_lines=file_lines, i_line=i_line,
                 key_counter=key_counter, line=line)): break
            if(len(line) == len(key_counter.keys)):
              try:
                if(key_counter.i_fobs is not None):
                  data_ = line[key_counter.i_fobs]
                  if(key_counter.i_sfobs is not None):
                    sigma_ = line[key_counter.i_sfobs]
                else:
                  if(key_counter.i_siobs is not None):
                    sigma_ = line[key_counter.i_siobs]
                  data_ = line[key_counter.i_iobs]
                if(key_counter.i_flag is not None):
                  flag_ = line[key_counter.i_flag]
                if(data_.count("*")>0 or data_.count("?")>0 or data_=="."):
                  continue
                if(sigma_ is not None):
                  if(sigma_.count("*")>0 or sigma_.count("?")>0 or sigma_=="."):
                    sigma_ = 1.0
              except:
                self.reset(message ="Cannot extract column data,#1.",line=line)
                break
              try:
                data_ = float(data_)
                if(data_ == 0.0): continue
                if(sigma_ is not None):
                  sigma_ = float(sigma_)
                  if(sigma_ < 0.0):
                    if(data_ == sigma_): continue
              except:
                if(self.is_break(file_lines=file_lines, i_line=i_line,
                   key_counter=key_counter, line=line)): break
              assert result_hkl.count(None) == 0
              assert data_ is not None
              if(result_hkl.count(0) != 3 and data_ != 0):
                if(max(max(result_hkl), abs(min(result_hkl))) > 10000):
                  self.reset(message ="Too big Miller index (> 10000).",line=line)
                  break
                self.indices.append(result_hkl)
                self.data.append(data_)
                if(flag_ is not None): self.flags.append(flag_)
                if(sigma_ is not None): self.sigmas.append(sigma_)
        else:
          if(self.data.size() > 0 and len(line) != len(key_counter.keys)):
            if(self.is_break(file_lines=file_lines, i_line=i_line,
               key_counter=key_counter, line=line)): break
    if(self.indices.size() != self.data.size()):
      self.reset(message = "self.indices.size() != self.data.size()")
    if(len(self.sigmas) > 0):
      if(self.indices.size() != self.sigmas.size()):
        self.reset(message = "self.indices.size() != self.sigmas.size()")
    else:
      self.sigmas = flex.double(self.data.size(), 1.0)
    if(self.indices.size() > 0 and self.flags.size() > 0):
      if(self.indices.size() != self.flags.size()):
        self.reset(message = "self.indices.size() != self.flags.size()")
      else:
        self.flags, cif_selection = get_array_of_r_free_flags(
          flags            = self.flags,
          crystal_symmetry = crystal_symmetry,
          indices          = self.indices,
          file_name        = file_name)
        if([self.flags, cif_selection].count(None) == 0):
          self.indices = self.indices.select(cif_selection)
          self.data = self.data.select(cif_selection)
          self.sigmas = self.sigmas.select(cif_selection)
          self.flags = self.flags.select(cif_selection)
    if(self.indices.size() == 0):
      print "No data extracted from input cif file."


  def access_hkl(self, line, key_counter):
    h, k, l = [None]*3
    try:
      h_ = line[key_counter.i_h]
      k_ = line[key_counter.i_k]
      l_ = line[key_counter.i_l]
      try_h = h_.replace("-","").strip().isdigit()
      try_k = k_.replace("-","").strip().isdigit()
      try_l = l_.replace("-","").strip().isdigit()
      if(not (try_h and try_k and try_l)): return h, k, l
      h = int(h_)
      k = int(k_)
      l = int(l_)
    except:
      return h, k, l
    return h, k, l

  def get_next_line(self, file_lines, i_line):
    next_line = None
    try:
      next_line = file_lines[i_line+1]
      if(len(next_line) == 0):
        try:
          next_line = file_lines[i_line+2]
        except: pass
    except: pass
    if(next_line is not None):
      next_line = next_line.strip()
      next_line = next_line.split()
    return next_line

  def is_break(self, file_lines, i_line, key_counter, line):
    next_line = self.get_next_line(file_lines=file_lines,i_line=i_line)
    if(next_line is not None):
      next_line_hkl = list(self.access_hkl(line=next_line,
        key_counter=key_counter))
      if(next_line_hkl.count(None) == 0):
        self.reset(message="Cannot extract column data,#2.",line=line)
        return True
      else:
        return True
    else: return True
    return False

  def reset(self, message, line=None):
    self.indices = flex.miller_index()
    self.data = flex.double()
    self.sigmas = flex.double()
    self.flags = flex.std_string()
    suffix = ""
    if(line is not None): suffix = line
    print message, self.file_name, suffix

def create_mtz_object(pre_miller_arrays,
                      key_counter,
                      crystal_symmetry,
                      file_name,
                      show_details_if_error):
  r_free_flags = None
  if(pre_miller_arrays.flags is not None):
    if(pre_miller_arrays.flags.size() > 0):
      r_free_flags = pre_miller_arrays.flags
  try:
    miller_set = miller.set(
      crystal_symmetry = crystal_symmetry,
      indices          = pre_miller_arrays.indices)
    miller_set = miller_set.auto_anomalous()
    miller_array = miller_set.array(
      data   = pre_miller_arrays.data,
      sigmas = pre_miller_arrays.sigmas)
    if(key_counter.i_fobs is not None):
      miller_array.set_observation_type_xray_amplitude()
      observation_type = "FOBS"
    else:
      miller_array.set_observation_type_xray_intensity()
      observation_type = "IOBS"
  except Exception, e:
    print "Cannot setup final miller array:", file_name, str(e)
    return None
  if(r_free_flags is not None):
    r_free_flags = miller_array.array(data = r_free_flags)
    assert r_free_flags.indices().all_eq(miller_array.indices())
  n_all = miller_array.indices().size()
  sel_unique = miller_array.unique_under_symmetry_selection()
  sel_dup = ~flex.bool(n_all, sel_unique)
  n_duplicate = sel_dup.count(True)
  n_uus = sel_unique.size()
  if(n_uus != n_all and n_duplicate > 1):
    print "Miller indices not unique under symmetry:", file_name, \
      "(%d redundant indices out of %d)" % (n_all-n_uus, n_all)
    if (show_details_if_error):
      miller_array.show_comprehensive_summary(prefix="  ")
      miller_array.map_to_asu().sort().show_array(prefix="  ")
    return None
  elif(n_duplicate == 1):
    miller_array = miller_array.select(sel_unique)
    if(r_free_flags is not None):
      r_free_flags = r_free_flags.select(sel_unique)
  mtz_dataset = miller_array.as_mtz_dataset(
    column_root_label = observation_type)
  if(r_free_flags is not None):
    mtz_dataset.add_miller_array(
      miller_array      = r_free_flags,
      column_root_label = "R-free-flags")
    assert r_free_flags.indices().all_eq(miller_array.indices())
  mtz_object = mtz_dataset.mtz_object()
  return mtz_object

class guess_observation_type(object):
  def __init__(self, file_name, pdb_raw_records, mtz_object, show_log):
    self.file_name = file_name
    self.mtz_object = mtz_object
    self.pdb_raw_records = pdb_raw_records
    xray_structure = self.get_xray_structure_from_pdb_file()
    if(xray_structure is not None):
      miller_arrays = self.mtz_object.as_miller_arrays()
      if(len(miller_arrays) == 1):
        f_obs = miller_arrays[0]
        r_free_flags = f_obs.array(data = flex.bool(f_obs.data().size(),False))
        r_free_flags_orig = None
      else:
        assert len(miller_arrays) == 2
        for miller_array in miller_arrays:
          if(miller_array.observation_type() is not None):
            f_obs = miller_array
          else:
            r_free_flags = miller_array
            r_free_flags_orig = r_free_flags.deep_copy()
            r_free_flags = r_free_flags.array(data = r_free_flags.data()==1)
      f_obs = f_obs.common_set(r_free_flags)
      r_free_flags = r_free_flags.common_set(f_obs)
      if(r_free_flags_orig is not None):
        f_obs = f_obs.common_set(r_free_flags_orig)
        r_free_flags_orig = r_free_flags_orig.common_set(f_obs)
      observation_type = f_obs.observation_type()
      if(str(observation_type) == "xray.amplitude"):
        observation_type = "F"
      elif(str(observation_type) == "xray.intensity"):
        observation_type = "I"
      else:
        assert observation_type is not None
      fmodel = self.get_fmodel_object(
        xray_structure = xray_structure,
        f_obs          = f_obs,
        r_free_flags   = r_free_flags)
      if(fmodel is not None):
        result, r_free_is_ok = self.get_observation_type(fmodel = fmodel,
          show_log = show_log)
        if([result, r_free_is_ok].count(None) == 0):
          observation_type = result
          # XXX do it optimally
          if(observation_type.startswith("F")):
            f_obs.set_observation_type_xray_amplitude()
          elif(observation_type.startswith("I")):
            f_obs.set_observation_type_xray_intensity()
          else:
            f_obs = f_obs.set_observation_type(observation_type = None)
          assert r_free_flags.indices().all_eq(f_obs.indices())
          if(r_free_flags_orig is not None):
            assert r_free_flags_orig.indices().all_eq(f_obs.indices())
          mtz_dataset=f_obs.as_mtz_dataset(column_root_label=observation_type)
          if(r_free_flags_orig is not None and r_free_is_ok):
            mtz_dataset.add_miller_array(
              miller_array      = r_free_flags_orig,
              column_root_label = "R-free-flags")
          self.mtz_object = mtz_dataset.mtz_object()

  def get_xray_structure_from_pdb_file(self):
    mon_lib_srv = monomer_library.server.server()
    ener_lib = monomer_library.server.ener_lib()
    xray_structure = None
    try:
      processed_pdb_file = monomer_library.pdb_interpretation.process(
        mon_lib_srv = mon_lib_srv,
        ener_lib    = ener_lib,
        raw_records = self.pdb_raw_records)
      xray_structure = processed_pdb_file.xray_structure()
    except Exception, e:
      print "INFO: Cannot extract xray structure:",self.file_name,str(e)
      return None
    if(xray_structure.scatterers().size()==0):
      print "INFO: Empty xray_structure:", self.file_name
      return None
    occ = xray_structure.scatterers().extract_occupancies()
    beq = xray_structure.extract_u_iso_or_u_equiv()
    sca = xray_structure.scatterers().extract_scattering_types()
    sel = (occ > 0.) & (occ < 50.) & (beq > 0.) & (beq < 999.) & (sca != "?")
    if(sel.size() > 0 and sel.count(True) > 0):
      xray_structure = xray_structure.select(selection = sel)
      xray_structure.tidy_us()
    else: return None
    return xray_structure

  def get_fmodel_object(self, xray_structure, f_obs, r_free_flags):
    fmodel = None
    selection = f_obs.d_spacings().data() > 0.25
    f_obs = f_obs.select(selection)
    r_free_flags = r_free_flags.select(selection)
    try:
      fmodel = mmtbx.f_model.manager(
        xray_structure   = xray_structure,
        r_free_flags     = r_free_flags,
        target_name      = "ls_wunit_k1",
        f_obs            = f_obs)
    except Exception, e:
      print "INFO: Cannot setup fmodel:", self.file_name, str(e)
      return None
    try:
      fmodel.update_xray_structure(
        xray_structure = xray_structure,
        update_f_calc  = True,
        update_f_mask  = True)
    except Exception, e:
      print "INFO: fmodel.update_xray_structure() failed:",\
        self.file_name,str(e)
      return None
    return fmodel

  def update_solvent_and_scale_helper(self, fmodel, params):
    try:
      fmodel.update_solvent_and_scale(params = params, verbose = -1)
      r_work = fmodel.r_work()*100
      r_free = fmodel.r_free()*100
    except Exception, e:
      print "INFO: fmodel.update_solvent_and_scale() failed:",\
        self.file_name,str(e)
      return None, None, None
    return fmodel, r_work, r_free

  def f_or_i_or_n(self,
                  fmodel,
                  observation_type,
                  params,
                  results):
    fmodel_dc = fmodel.deep_copy()
    observation_type_dc = observation_type
    f_obs = fmodel_dc.f_obs
    f_obs = f_obs.set_observation_type(observation_type = None)
    f_obs = f_obs.f_sq_as_f()
    fmodel_dc.update(f_obs = f_obs)
    fmodel_, r_work_, r_free_ = self.update_solvent_and_scale_helper(
      fmodel = fmodel_dc, params = params)
    if(fmodel_ is not None):
      observation_type_dc = "I"+observation_type_dc[1:]
      results.append([r_work_, observation_type_dc, r_free_])
    #
    fmodel_dc = fmodel.deep_copy()
    observation_type_dc = observation_type
    f_obs = fmodel_dc.f_obs
    f_obs = f_obs.set_observation_type(observation_type = None)
    f_obs = f_obs.f_as_f_sq()
    fmodel_dc.update(f_obs = f_obs)
    fmodel_, r_work_, r_free_ = self.update_solvent_and_scale_helper(
      fmodel = fmodel_dc, params = params)
    if(fmodel_ is not None):
      observation_type_dc = "F"+observation_type_dc[1:]
      results.append([r_work_, observation_type_dc, r_free_])
    return results

  def get_observation_type(self, fmodel, show_log):
    params = bss.master_params.extract()
    params.k_sol_max = 0.6
    params.k_sol_min = 0.0
    params.b_sol_max = 500.0
    params.b_sol_min = 0.0
    neutron_flag = True
    observation_type = str(fmodel.f_obs.observation_type())
    assert observation_type in ["xray.amplitude","xray.intensity"]
    if(observation_type == "xray.amplitude"): observation_type = "FOBS_X"
    if(observation_type == "xray.intensity"): observation_type = "IOBS_X"
    observation_type_from_file = observation_type
    d_max, d_min = fmodel.f_obs.d_max_min()
    selection = fmodel.f_obs.all_selection()
    if(d_min < 3.0):
      selection &= fmodel.f_obs.d_spacings().data() <= 6.0
    selection &= fmodel.f_obs.d_spacings().data() >= 2.0
    sigmas = fmodel.f_obs.sigmas()
    if(abs(flex.mean(sigmas) - 1.0) > 0.001):
      selection &= (fmodel.f_obs.data() > fmodel.f_obs.sigmas()*2.0)
    selection &= (fmodel.f_obs.d_star_sq().data() > 0)
    ssize = selection.size()
    if(ssize > 0 and selection.count(True) > 100):
      fmodel = fmodel.select(selection = selection)
    try:
      fmodel = fmodel.remove_outliers()
    except Exception, e:
      print "INFO: fmodel.remove_outliers() failed:",self.file_name,str(e)
    results = []
    #
    fmodel_dc = fmodel.deep_copy()
    fmodel_dc, r_work, r_free = self.update_solvent_and_scale_helper(
      fmodel = fmodel_dc, params = params)
    if(fmodel_dc is not None):
      results.append([r_work, observation_type, r_free])
    if(r_work > 30.0):
      results = self.f_or_i_or_n(fmodel = fmodel_dc, params = params,
        observation_type = observation_type, results = results)
      #
      fmodel_dc = fmodel.deep_copy()
      # XXX push to xray_structure ?
      neutron_scattering_dict = {}
      reg = fmodel_dc.xray_structure.scattering_type_registry()
      for scattering_type in reg.type_index_pairs_as_dict().keys():
        scattering_info = neutron_news_1992_table(scattering_type, 1)
        b = scattering_info.bound_coh_scatt_length()
        if(b.imag != 0.0):
          neutron_flag = False
          break
        neutron_scattering_dict[scattering_type] = \
          eltbx.xray_scattering.gaussian(b.real)
      if(neutron_flag):
        fmodel_dc.xray_structure.scattering_type_registry(
          custom_dict = neutron_scattering_dict)
        fmodel_dc.update_xray_structure(update_f_calc=True, update_f_mask=True)
        fmodel_dc, r_work, r_free = self.update_solvent_and_scale_helper(
          fmodel = fmodel_dc, params = params)
        observation_type = observation_type[:-2] + "_N"
        results.append([r_work, observation_type, r_free])
        results = self.f_or_i_or_n(fmodel = fmodel_dc, params = params,
          observation_type = observation_type, results = results)
    try:
      result_best = observation_type
      r_free_is_ok = False
      if show_log: print "Scores:"
      if show_log:
        print "  observation type from input file:", observation_type_from_file
      if(len(results) == 1):
        if(results[0][0] < 30.):
          observation_type = "F"+observation_type[1:]
          r_free_is_ok = (results[0][0] < results[0][2]) and \
            (abs(results[0][0] - results[0][2]) > 1.0)
      elif(len(results) > 1):
        r = results[0][0]
        for res in results:
          if show_log: print "    if it is %s: r_work= %6.2f"%(res[1], res[0])
          if(res[0] <= r):
            result_best = res
            r = res[0]
        delta = 100.
        for res in results:
          if(res[0] > result_best[0] and abs(res[0]-result_best[0]) < delta):
            delta = abs(res[0]-result_best[0])
        if show_log:
          print "  selected %s with r_work %6.2f"%(result_best[1],result_best[0])
        r_free_is_ok = (result_best[0] < result_best[2]) and \
          (abs(result_best[0] - result_best[2]) > 1.0)
        if((delta < 3. and result_best[0] > 30.) or result_best[0] > 40.):
          observation_type = "OBS"
        else: observation_type = result_best[1]
        if(observation_type != "OBS"):
          if(neutron_flag and results[3][0] == result_best[0]):
            observation_type = "F"+observation_type[1:]
          if(results[0][0] == result_best[0]):
            observation_type = "F"+observation_type[1:]
      if show_log: print "  final observation type:", observation_type
      if(len(results) == 1):
        if show_log: print "  final r_work= %6.2f"%results[0][0]
    except Exception, e:
      print "INFO: Cannot score results to select observation type:", \
        self.file_name,str(e)
      return None, None
    return observation_type, r_free_is_ok

def run(args, command_name = "phenix.cif_as_mtz"):
  if (len(args) == 0): args = ["--help"]
  try:
    command_line = (iotbx_option_parser(
      usage="%s [reflection_cif_file] [options]" % command_name,
      description='Example: %s r1o9ksf.ent --symmetry=pdb1o9k.ent'%command_name)
      .enable_symmetry_comprehensive()
      .option(None, "--output_file_name",
        action="store",
        default=False,
        type="string",
        help="Output mtz file name.")
      .option(None, "--use_model",
        action="store",
        default=False,
        type="string",
        help="Use PDB model to make better guess about reflection data type.")
      .option(None, "--wavelength_id",
        action="store",
        default=None,
        type="int",
        help="Extract data set with given wavelength_id.")
      .option(None, "--crystal_id",
        action="store",
        default=None,
        type="int",
        help="Extract data set with given crystal_id.")
      .option("--show_details_if_error",
          action="store_true",
          help="Show data details for some errors.")
      .option("--show_log",
          action="store_true",
          help="Show some output.")
    ).process(args=args)
  except Exception, e:
    if(str(e) != "0"): print str(e)
    sys.exit(0)
  crystal_symmetry = command_line.symmetry
  if(command_line.symmetry.unit_cell() is None or
     command_line.symmetry.space_group_info() is None):
    if(command_line.options.use_model):
      crystal_symmetry = crystal_symmetry_from_pdb.extract_from(
         file_name=command_line.options.use_model)
  if(crystal_symmetry.unit_cell() is None or
     crystal_symmetry.space_group_info() is None):
    raise Sorry(
      "Crystal symmetry is not defined. Please use the --symmetry option.\n"
      "Type %s without arguments to see more options."%command_name)
  if(len(command_line.args) > 1):
    print "%d arguments are given from the command line:"% \
      len(command_line.args), command_line.args
    raise Sorry("Please specify one reflection cif file.")
  file_name = command_line.args[0]
  if(not os.path.isfile(file_name)):
    raise Sorry("File is not found: %s"%file_name)
  file_lines = smart_open.for_reading(file_name=file_name).read().splitlines()
  mtz_object = extract(
    file_name             = file_name,
    file_lines            = file_lines,
    crystal_symmetry      = crystal_symmetry,
    wavelength_id         = command_line.options.wavelength_id,
    crystal_id            = command_line.options.crystal_id,
    show_details_if_error = command_line.options.show_details_if_error)
  if(mtz_object is not None):
    pdb_file_name = command_line.options.use_model
    if(pdb_file_name):
      pdb_raw_records = smart_open.for_reading(
        file_name=pdb_file_name).read().splitlines()
      mtz_object = guess_observation_type(
        file_name       = pdb_file_name,
        pdb_raw_records = pdb_raw_records,
        mtz_object      = mtz_object,
        show_log        = command_line.options.show_log).mtz_object
    if(command_line.options.output_file_name):
      output_file_name = command_line.options.output_file_name
    else:
      basename = os.path.basename(file_name)
      if(basename[-4:-3] == "."): output_file_name = basename[:-4]+".mtz"
      elif(basename[-5:-4] == "."): output_file_name = basename[:-5]+".mtz"
      elif(basename.endswith(".ent.gz")): output_file_name=basename[:-7]+".mtz"
      else: output_file_name = basename+".mtz"
    mtz_object.write(file_name = output_file_name)

def extract(file_name, file_lines, crystal_symmetry, wavelength_id, crystal_id,
            show_details_if_error):
  keys = extract_keys(file_name = file_name, file_lines = file_lines)
  if(len(keys) == 0): return None
  key_counter = count_keys(keys = keys, file_name = file_name)
  if(wavelength_id is not None):
    if("wavelength_id" not in keys or key_counter.i_wavelength_id is None):
      raise Sorry("Input cif file does not contain wavelength_id key.")
  if(crystal_id is not None):
    if("crystal_id" not in keys or key_counter.i_crystal_id is None):
      raise Sorry("Input cif file does not contain crystal_id key.")
  if(key_counter.i_h is None): return None
  pre_miller_arrays = extract_data(
    key_counter       = key_counter,
    file_name         = file_name,
    file_lines        = file_lines,
    wavelength_id     = wavelength_id,
    crystal_id        = crystal_id,
    crystal_symmetry  = crystal_symmetry)
  if(pre_miller_arrays.indices.size() == 0): return None
  mtz_object = create_mtz_object(
    pre_miller_arrays = pre_miller_arrays,
    key_counter       = key_counter,
    crystal_symmetry  = crystal_symmetry,
    file_name         = file_name,
    show_details_if_error = show_details_if_error)
  return mtz_object

if(__name__ == "__main__"):
   run(sys.argv[1:])
