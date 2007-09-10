# LIBTBX_SET_DISPATCHER_NAME phenix.cif_as_mtz

import sys, os, time
from cctbx.array_family import flex
from cctbx import miller
from libtbx import easy_run
from cctbx import crystal
from iotbx.option_parser import iotbx_option_parser
from libtbx.utils import Sorry
from iotbx import reflection_file_utils

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

ccp4_range = range(-99,-1)+ range(2,150)
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
    n_h = 0
    n_k = 0
    n_l = 0
    n_fobs = 0
    n_iobs = 0
    n_flags = 0
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
                     crystal_symmetry):
    self.indices = flex.miller_index()
    self.data = flex.double()
    self.sigmas = flex.double()
    self.flags = flex.std_string()
    self.file_name = file_name
    start_counter = 0
    start_flag = False
    for line in file_lines:
      line_orig = line
      h_,k_,l_,data_,sigma_,flag_ = [None]*6
      line = line.strip()
      line = line.split()
      if(len(key_counter.keys) != start_counter or not start_flag):
        if(len(line) == 1):
          if(line[0] == "loop_"): start_flag = True
          if(start_flag and (line[0].replace("_refln.","") in key_counter.keys)
             and line[0] != "loop_"):
            start_counter += 1
      if(len(key_counter.keys) == start_counter and start_flag):
        if(self.indices.size() > 0 and len(line) != len(key_counter.keys)):
          break
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
            if([data_,sigma_,flag_].count("?") > 0): continue
          except:
            self.reset(message = "Cannot extract column data,#1.",line=line)
            break
          try:
            h_ = line[key_counter.i_h]
            k_ = line[key_counter.i_k]
            l_ = line[key_counter.i_l]
            try_h = h_.replace("-","").strip().isdigit()
            try_k = k_.replace("-","").strip().isdigit()
            try_l = l_.replace("-","").strip().isdigit()
            if(not (try_h and try_k and try_l)):
              if(line[0] not in ["#END","#","#;"] and
                 line_orig.count("This file should contain")==0):
                self.reset(message="h or k or k or all: NOT digit(s).",line=line)
              break
            h_ = int(h_)
            k_ = int(k_)
            l_ = int(l_)
            data_ = float(data_)
            if(data_ == 0.0): continue
            if(sigma_ is not None): sigma_ = float(sigma_)
          except:
            self.reset(message = "Cannot extract column data,#2.",line=line)
            break
          assert [h_, k_, l_].count(None) == 0
          assert data_ is not None
          if([h_,k_,l_].count(0) != 3 and data_ != 0):
            self.indices.append([h_, k_, l_])
            self.data.append(data_)
            if(flag_ is not None): self.flags.append(flag_)
            if(sigma_ is not None): self.sigmas.append(sigma_)
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
                      file_name):
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
  try:
    mtz_dataset = miller_array.as_mtz_dataset(
      column_root_label = observation_type)
  except Exception, e:
    print "Cannot do miller_array.as_mtz_dataset():", str(e), file_name
    return None
  if(pre_miller_arrays.flags is not None):
    if(pre_miller_arrays.flags.size() > 0):
      mtz_dataset.add_miller_array(
        miller_array      = miller_array.array(data = pre_miller_arrays.flags),
        column_root_label = "R-free-flags")
  mtz_object = mtz_dataset.mtz_object()
  return mtz_object

def run(args, command_name = "phenix.cif_as_mtz"):
  if (len(args) == 0): args = ["--help"]
  try:
    command_line = (iotbx_option_parser(
      usage="%s [reflection_cif_file] [options]" % command_name,
      description='Example: %s r1o9ksf.ent --symmetry=pdb1o9k.ent'%command_name)
      .enable_show_defaults()
      .enable_symmetry_comprehensive()
      .option("--remove_input_file",
          action="store_true",
          help="Remove input CIF file (very dangerous option).")
      .option(None, "--output_file_name",
        action="store",
        default=False,
        type="string",
        help="Output mtz file name.")
    ).process(args=args)
  except Exception, e:
    if(str(e) != "0"): print str(e)
    sys.exit(0)
  crystal_symmetry = command_line.symmetry
  if(command_line.symmetry.unit_cell() is None or
     command_line.symmetry.space_group_info() is None):
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
  ifo = open(file_name, "r")
  file_lines = ifo.readlines()
  ifo.close()
  mtz_object = extract(
    file_name        = file_name,
    file_lines       = file_lines,
    crystal_symmetry = crystal_symmetry)
  if(mtz_object is not None):
    if(command_line.options.output_file_name):
      output_file_name = command_line.options.output_file_name
    else:
      if(file_name[-4:-3] == "."): output_file_name = file_name[:-4]+".mtz"
      elif(file_name[-5:-4] == "."): output_file_name = file_name[:-5]+".mtz"
      else: output_file_name = file_name+".mtz"
    mtz_object.write(file_name = output_file_name)
  if(command_line.options.remove_input_file):
    os.remove(file_name)

def extract(file_name, file_lines, crystal_symmetry):
  keys = extract_keys(file_name = file_name, file_lines = file_lines)
  if(len(keys) == 0): return None
  key_counter = count_keys(keys = keys, file_name = file_name)
  if(key_counter.i_h is None): return None
  pre_miller_arrays = extract_data(
    key_counter       = key_counter,
    file_name         = file_name,
    file_lines        = file_lines,
    crystal_symmetry  = crystal_symmetry)
  if(pre_miller_arrays.indices.size() == 0): return None
  mtz_object = create_mtz_object(
    pre_miller_arrays = pre_miller_arrays,
    key_counter       = key_counter,
    crystal_symmetry  = crystal_symmetry,
    file_name         = file_name)
  return mtz_object

if(__name__ == "__main__"):
   run(sys.argv[1:])
