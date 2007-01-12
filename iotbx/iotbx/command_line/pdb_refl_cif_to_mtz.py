import sys, os, time
from cctbx.array_family import flex
from iotbx import crystal_symmetry_from_any
from cctbx import miller
from iotbx import reflection_file_utils
import time, sys, os
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from iotbx import reflection_file_reader
import mmtbx.f_model
import mmtbx.bulk_solvent.bulk_solvent_and_scaling as bss
from mmtbx.scaling import outlier_rejection


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

ccp4_range = range(-99,-1)+ range(2,100)
tw_flags_num = ["1","1.","1.0","1.00","1.000","1.0000","-1","-1.","-1.0","-1.00",
            "-1.000","0","0.0","0.","0.00","0.000","0.0000"]+ \
    [str(i) for i in ccp4_range] + ccp4_range + \
    [str("%.1f"%i) for i in ccp4_range]+[str("%.2f"%i) for i in ccp4_range]+\
    [str("%.3f"%i) for i in ccp4_range]+[str("%.1f"%i)[:-1] for i in ccp4_range]

cif_flags = ["x","o","f","-",">","<","h","l","O","F"]

tw_flags_all = tw_flags_num + cif_flags

###############################################################################
mon_lib_srv = monomer_library.server.server()
ener_lib    = monomer_library.server.ener_lib()
line = "PDBid  d_min  d_max  k_sol  b_sol  r_work  r_free       B11     B22     B33     B12     B13     B23   SG    model     data  flags_%"

def get_xray_structure_from_pdb_file(pdb, out):
  try:
    processed_pdb_file = monomer_library.pdb_interpretation.process(
                                                     mon_lib_srv = mon_lib_srv,
                                                     ener_lib    = ener_lib,
                                                     file_name   = pdb)
    xray_structure = processed_pdb_file.xray_structure()
  except KeyboardInterrupt: raise
  except:
    xray_structure = None
    print >> out, "ERROR_(cannot extract xray_structure): %s"%pdb[-11:]
  try:
    if(xray_structure is not None and xray_structure.scatterers().size()>0):
       occ = xray_structure.scatterers().extract_occupancies()
       beq = xray_structure.extract_u_iso_or_u_equiv()
       sel = (occ > 0.0) | (occ < 5.0) | (beq > 0.0) | (beq < 500.0)
       xray_structure = xray_structure.select(selection = sel)
  except KeyboardInterrupt: raise
  except:
    xray_structure = None
    print >> out, "ERROR_(cannot extract xray_structure): %s"%pdb[-11:]
  if(xray_structure is not None and xray_structure.scatterers().size()==0):
    xray_structure = None
    print >> out, "ERROR_(cannot extract xray_structure): %s"%pdb[-11:]
  return xray_structure

def get_fobs_and_rfree_flags_from_mtz(hkl, pdb, out):
  try:
    f_obs, r_free_flags = None, None
    refl = reflection_file_reader.any_reflection_file(file_name = hkl)
    refl_arrays = refl.as_miller_arrays()
    for a in refl_arrays:
        if(str(a.info()).count("FOBS")):
           f_obs = a
           assert str(f_obs.observation_type()) == "xray.amplitude"
        if(str(a.info()).count("FLAGS")):
           r_free_flags = a
        if(str(a.info()).count("IOBS")):
           assert str(a.observation_type()) == "xray.intensity"
           f_obs = a.f_sq_as_f()
    r_free_flags = r_free_flags.array(data = (r_free_flags.data() == 1))
    b, a = r_free_flags.common_sets(f_obs)
    if(a.indices().size() != f_obs.indices().size()):
       assert a.indices().size() == f_obs.indices().size()
    else:
       r_free_flags = b
       assert a.indices().all_eq(f_obs.indices())
       assert str(f_obs.observation_type()) == "xray.amplitude"
  except KeyboardInterrupt: raise
  except:
    f_obs, r_free_flags = None, None
    print >> out, "ERROR_(inconsistency in f_obs or/and r_free_flags): %s %s"%(hkl[:-4],pdb)
  return f_obs, r_free_flags

def get_output_file_object():
  t = time.ctime().split()
  time_stamp = t[4]+"_"+t[1].upper()+"_"+t[2]+"_"+t[3][:-3].replace(":","h")
  return open("pdb_ksol_bsol_"+time_stamp, "w"), \
         open("errors_pdb_ksol_bsol_"+time_stamp, "w")

def get_fmodel_object(xray_structure, f_obs, r_free_flags, pdb, hkl, out):
  try:
    fmodel = mmtbx.f_model.manager(xray_structure   = xray_structure,
                                   sf_algorithm     = "fft",
                                   sf_cos_sin_table = True,
                                   r_free_flags     = r_free_flags,
                                   target_name      = "ls_wunit_k1",
                                   f_obs            = f_obs)
    fmodel.update_xray_structure(xray_structure = xray_structure,
                                 update_f_calc  = True,
                                 update_f_mask  = True)
  except Exception, e:
    print >> out, "ERROR_(problem in setting up fmodel): %s"%pdb[-11:], hkl[:-4], str(e)
    fmodel = None
  return fmodel

def update_solvent_and_scale_helper(fmodel, params, low_res = 6.0,
                              sigma_cutoff = 2.0, out=None,pdb=None, hkl=None):
  proceed = True
  r_work, r_free = None,None
  f_obs = fmodel.f_obs
  selection = f_obs.all_selection()
  selection &= f_obs.d_spacings().data() <= low_res
  selection &= (f_obs.data() > f_obs.sigmas()*sigma_cutoff)
  selection &= (f_obs.d_star_sq().data() > 0)
  try:
    assert selection.count(True) > 0
  except Exception, e:
    print >> out, "ERROR_(No_rflections_left after sigma and resolutions cutoff:", pdb[-11:], hkl[:-4]
    proceed = False
  if(proceed):
     try:
       fmodel = fmodel.select(selection = selection)
       fmodel.update(k_sol=0, b_sol=0, b_cart=[0,0,0,0,0,0])
       fmodel.update_solvent_and_scale(params = params, verbose = -1)
       r_work = fmodel.r_work()*100
       r_free = fmodel.r_free()*100
     except Exception, e:
       print >> out, "ERROR_(update_solvent_and_scale2):", pdb[-11:], hkl[:-4]
  return r_work, r_free

def update_solvent_and_scale(fmodel, pdb, hkl, out):
  try:
    status = None
    params = bss.solvent_and_scale_params()
    params.b_sol_min=0.0
    params.b_sol_max=200.0
    params.k_sol_min=0.0
    params.k_sol_max=1.5
    fmodel.update_solvent_and_scale(params = params, verbose = -1)
    r_work = fmodel.r_work()*100
    r_free = fmodel.r_free()*100
    rwrf_delta = 1.5
    if(fmodel.f_obs.d_min() < 1.2): rwrf_delta = 0.5
    if(r_work > r_free or abs(r_work - r_free) <= rwrf_delta or r_work > 45.0):
       status = "bad"
    else:
       status = "good"
  except KeyboardInterrupt: raise
  except:
    print >> out, "ERROR_(update_solvent_and_scale1): %s"%pdb[-11:], hkl[:-4]
    status = None
    fmodel = None
  convert_to_intensities = False
  convert_to_amplitudes  = False
  if([status, fmodel].count(None) == 0 and status == "bad"):
     fmodel_dc = fmodel.deep_copy()
     rw, rf = update_solvent_and_scale_helper(fmodel = fmodel_dc,
                                      params = params, out=out,pdb=pdb,hkl=hkl)
     if([rw, rf].count(None) > 0): return None, None
     check = (rw>rf or abs(rw-rf) <= rwrf_delta or rw > 45.)
     if(check): status = "bad"
     else:      status = "good"
  if([status, fmodel].count(None) == 0 and status == "bad"):
     fmodel_dc = fmodel.deep_copy()
     f_obs = fmodel_dc.f_obs
     f_obs = f_obs.set_observation_type(observation_type = None)
     f_obs = f_obs.f_sq_as_f()
     fmodel_dc.update(f_obs = f_obs)
     rw, rf = update_solvent_and_scale_helper(fmodel = fmodel_dc,
                                      params = params, out=out,pdb=pdb,hkl=hkl)
     if([rw, rf].count(None) > 0): return None, None
     check = (rw>rf or abs(rw-rf) <= rwrf_delta or rw > 45. or abs(r_work-rw) <= 5.)
     if(check): status = "bad"
     else:
        status = "good"
        convert_to_amplitudes = True
  if([status, fmodel].count(None) == 0 and status == "bad"):
     fmodel_dc = fmodel.deep_copy()
     f_obs = fmodel_dc.f_obs
     f_obs = f_obs.set_observation_type(observation_type = None)
     f_obs = f_obs.f_as_f_sq()
     fmodel_dc.update(f_obs = f_obs)
     rw, rf = update_solvent_and_scale_helper(fmodel = fmodel_dc,
                                      params = params, out=out,pdb=pdb,hkl=hkl)
     if([rw, rf].count(None) > 0): return None, None
     check = (rw>rf or abs(rw-rf) <= rwrf_delta or rw > 45. or abs(r_work-rw) <= 5.)
     if(check): status = "bad"
     else:
        status = "good"
        convert_to_intensities = True
  if(status == "good" and
            [convert_to_intensities, convert_to_amplitudes].count(True) > 0):
     fmodel.f_obs.set_observation_type(observation_type = None)
     if(convert_to_intensities):
        fmodel.update(k_sol=0, b_sol=0, b_cart=[0,0,0,0,0,0],
                                            f_obs = fmodel.f_obs.f_as_f_sq())
     if(convert_to_amplitudes):
        fmodel.update(k_sol=0, b_sol=0, b_cart=[0,0,0,0,0,0],
                                            f_obs = fmodel.f_obs.f_sq_as_f())
     fmodel.f_obs.set_observation_type_xray_amplitude()
     fmodel.update_solvent_and_scale(params = params, verbose = -1)
  if(status == "bad"):
     print >> out, "ERROR_(status=bad): %s"%pdb[-11:], hkl[:-4]
     status = None
     fmodel = None
  return status, fmodel

def one_run(pdb, hkl, error_out):
  result = None
  fmodel = None
  format="%5s%5s%7.2f%7.2f%7.2f%7.2f%8.2f%8.2f  %8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%5d%9d%9d%9.2f"
  xray_structure = get_xray_structure_from_pdb_file(pdb = pdb, out = error_out)
  if(xray_structure is not None and xray_structure.scatterers().size()>0):
     f_obs, r_free_flags = get_fobs_and_rfree_flags_from_mtz(
                                         hkl = hkl, pdb = pdb, out = error_out)
     if([f_obs, r_free_flags].count(None) == 0):
        fmodel = get_fmodel_object(xray_structure = xray_structure,
                                   f_obs          = f_obs,
                                   r_free_flags   = r_free_flags,
                                   pdb            = pdb,
                                   hkl            = hkl,
                                   out = error_out)
        if(fmodel is not None):
           status, fmodel = update_solvent_and_scale(
                                 fmodel = fmodel, pdb = pdb, hkl = hkl, out = error_out)
           if([status, fmodel].count(None) == 0 and status == "good"):
              d_max, d_min = fmodel.f_obs.d_max_min()
              b0,b1,b2,b3,b4,b5 = fmodel.b_cart()
              sg = xray_structure.space_group().type().number()
              hkl_size = fmodel.f_obs.indices().size()
              pdb_size = fmodel.xray_structure.scatterers().size()
              flags_pc = fmodel.r_free_flags.data().count(True)*100./fmodel.r_free_flags.data().size()
              result = format%(hkl[:-4],pdb[-8:-4],d_min,d_max,fmodel.k_sol(),
                fmodel.b_sol(), fmodel.r_work()*100, fmodel.r_free()*100,b0,
                b1,b2,b3,b4,b5,sg, pdb_size, hkl_size, flags_pc)
  return result, fmodel, format

###############################################################################

def get_array_of_r_free_flags(proceed, flags, cs, indices, pdb, err):
  pdb = pdb[-11:]
  flag_values = []
  guess_cif = 0
  cif_selection = None
  result_flags = flex.int()
  for item in flags:
      if(item not in tw_flags_all):
         print >> err, "Unknown r-free flag symbol: %s"%str(item), pdb, tw_flags_all
         proceed = False
         break
      if(item not in flag_values):
         flag_values.append(item)
  if(proceed):
     for item in flag_values:
         if(item in cif_flags): guess_cif += 1
     if(guess_cif > 0):
        cif_selection = flex.bool()
     for item in flags:
         item = item.strip()
         if(guess_cif == 0):
            if(str(item) in ["o","O"]): item = 0
            result_flags.append(int(float(item)))
         if(guess_cif > 0):
            if(item in ["f","F"]):
               result_flags.append(1)
               cif_selection.append(True)
            elif(item in ["o","O"]):
               result_flags.append(0)
               cif_selection.append(True)
            elif(item in ["-","x"]):
               result_flags.append(0)
               cif_selection.append(False)
            elif(item in tw_flags_num):
               result_flags.append(int(float(item)))
               cif_selection.append(True)
            else:
               result_flags.append(0)
               cif_selection.append(True)
     try:
       flags_mi = miller.set(
                   crystal_symmetry = cs,
                   indices          = indices).array(data = result_flags)
     except Exception, e:
       proceed = False
       print >> err, "Cannot_setup miller array with r-free flags: ",pdb,str(e)
     if(proceed):
        try:
          test_flag_value = reflection_file_utils.get_r_free_flags_scores(
                                    miller_arrays   = [flags_mi],
                                    test_flag_value = None).test_flag_values[0]
          if(test_flag_value is not None):
             result_flags = (result_flags == test_flag_value)
          else: proceed = False
        except Exception, e:
          proceed = False
          print >> err, "Cannot_score r-free flags1: ", pdb,str(e)
        if(test_flag_value is None):
           print >> err, "Cannot_score r-free flags2: ", pdb
  return result_flags, cif_selection, proceed

def compose_and_write_mtz(indices, data, sigmas, flags, cif_selection, pdb,
                                           err, proceed, cs, observation_type):
  file_name = None
  pdb = pdb[-11:]
  try:
    assert indices.size() == flags.size()
    assert sigmas.size() == data.size()
    assert data.size() == flags.size()
    assert flags.count(True) < flags.count(False)
    if(cif_selection is not None):
       indices = indices.select(cif_selection)
       sigmas  = sigmas.select(cif_selection)
       data    = data.select(cif_selection)
       flags   = flags.select(cif_selection)
    assert indices.size() == flags.size()
    assert sigmas.size() == data.size()
    assert data.size() == flags.size()
    assert flags.count(True) < flags.count(False)
    assert observation_type in ["FOBS", "IOBS"]
  except Exception, e:
    print >> err, pdb, str(e)
    proceed = False
  if(proceed):
     try:
       miller_array = miller.set(
                        crystal_symmetry = cs,
                        indices          = indices,
                        anomalous_flag   = False).array(data   = data,
                                                        sigmas = sigmas)
       if(observation_type == "FOBS"):
          miller_array.set_observation_type_xray_amplitude()
       elif(observation_type == "IOBS"):
          miller_array.set_observation_type_xray_intensity()
     except Exception, e:
       proceed = False
       print >> err, "Cannot_setup miller array: ", pdb, str(e)
     if(proceed):
        try:
          mtz_dataset = miller_array.as_mtz_dataset(
                                          column_root_label = observation_type)
        except Exception, e:
          proceed = False
          print >> err, "Cannot_write MTZ file= ", str(e), pdb
        if(proceed):
           mtz_dataset.add_miller_array(
                          miller_array      = miller_array.array(data = flags),
                          column_root_label = "FLAGS")
           mtz_object = mtz_dataset.mtz_object()
           file_name = pdb[3:-4]+".mtz"
           mtz_object.write(file_name = file_name)
  return proceed, file_name


def run(args):
  pdb = args[0]
  hkl = args[1]
  err = open("_error_"+pdb[-8:-4],"w")
  proceed = True
  data    = flex.double()
  sigmas  = flex.double()
  indices = flex.miller_index()
  flags_  = []
  flags = flex.bool()
  try:
    cs = crystal_symmetry_from_any.extract_from(pdb)
  except Exception, e:
    proceed = False
    print >> err, "Cannot_extract_crystal symmetry: ", str(e), pdb
  if(cs is None):
     proceed = False
     print >> err, "Cannot_extract_crystal symmetry: ", pdb
  #
  # getting keys
  #
  if(hkl.count(".Z") == 1 and proceed):
     os.system("cp "+ hkl + " .")
     hkl = hkl[-13:]
     #os.system("$HOME/uncompress " + hkl) # to work on BSGC-2
     os.system("gunzip " + hkl)
     print hkl
     ifo = open(hkl[:-2], "r")
     lines = ifo.readlines()
     ifo.close()
     # getting keys
     lock_0 = False
     lock_1 = False
     keys = []
     number_of_loop_records = 0
     for i_seq, st in enumerate(lines):
        st = st.strip()
        if(st.startswith("loop_")):
           if(i_seq < len(lines)):
              if(lines[i_seq+1].strip().startswith("_refln.")):
                 lock_0 = True
        if(lock_0 and st.startswith("_refln.")):
           lock_1 = True
        if(lock_0 and lock_1 and not st.startswith("_refln.")):
           lock_0 = False
           lock_1 = False
        if(lock_0 and lock_1):
           keys.append(st)
        if(len(keys) > 0 and st.startswith("loop_")):
           print >> err, "Multiple keys: ", hkl, keys
           keys = []
           proceed = False
           break
     if(proceed and len(keys) <= 4):
        print >> err, "No keys: ", keys, hkl
        proceed = False
     #
     # we have keys now
     #
     if(proceed):
        h_i  = None ; n_h = 0
        k_i  = None ; n_k = 0
        l_i  = None ; n_l = 0
        f_i  = None ; n_f = 0
        i_i  = None ; n_i = 0
        sf_i = None
        si_i = None
        t_i  = None ; n_t = 0
        keys_updated = []
        for i_seq, key in enumerate(keys):
            key = key.strip()
            key = key.replace("_refln.","")
            keys_updated.append(key)
        for i_seq, key in enumerate(keys_updated):
            if(key == "index_h"):
               h_i = i_seq
               n_h += 1
            if(key == "index_k"):
               k_i = i_seq
               n_k += 1
            if(key == "index_l"):
               l_i = i_seq
               n_l += 1
            if(key in keys_f_obs):
               f_i = i_seq
               n_f += 1
            if(key in keys_i_obs):
               i_i = i_seq
               n_i += 1
            if(key in keys_sf_obs): sf_i = i_seq
            if(key in keys_si_obs): si_i = i_seq
            if(key in keys_status):
               t_i = i_seq
               n_t += 1
            if(key not in possible_keys):
               proceed = False
               print >> err, "Unknown_key= ", key,  keys_updated, hkl
               break
        if(n_t > 1 and proceed):
           print >> err, "Multiple CV sets: ", keys_updated, hkl
           proceed = False
        if(n_t == 0 and proceed):
           proceed = False
           print >> err, "No CV sets: ", keys_updated, hkl, n_t
        if([f_i,i_i].count(None) == 2 or n_f > 1 or n_i > 1 and proceed):
           print >> err, "Multiple/no data sets: ", keys_updated,hkl,n_f,n_i
           proceed = False
        if([h_i,k_i,l_i].count(None) > 0 or n_h > 1 or n_k > 1 or n_l > 1 and proceed):
           print >> err, "Multiple/no index sets: ",\
                  keys_updated,hkl,n_h,n_k,n_l,"=",h_i,k_i,l_i
           proceed = False
     #
     # Hopefully, we know the numbers of colums with data. Try get it now.
     # Here we know the data is unique.
        if(proceed):
           for st in lines:
              if(st == ''): break
              if(proceed):
                 st = st.split()
                 if(indices.size() > 0 and len(st) != len(keys_updated)):
                    break
                 if(len(st) == len(keys_updated)):
                    try:
                      try_h = st[h_i].replace("-","").strip().isdigit()
                      try_k = st[k_i].replace("-","").strip().isdigit()
                      try_l = st[l_i].replace("-","").strip().isdigit()
                      is_digits = try_h and try_k and try_l
                      h_ = int(st[h_i])
                      k_ = int(st[k_i])
                      l_ = int(st[l_i])
                      if(f_i is not None): data_ = float(st[f_i])
                      else:                data_ = float(st[i_i])
                    except:
                      is_digits = False
                    if(is_digits and [h_,k_,l_].count(0)<3):
                       if(abs(h_)>10000 or abs(k_)>10000 or abs(l_)>10000):
                          proceed = False
                       indices.append([h_, k_, l_])
                       data.append(data_)
                       if(sf_i is not None):
                          assert f_i is not None
                          try:
                            sigmas.append(float(st[sf_i]))
                          except:
                            proceed = False
                       elif(si_i is not None):
                          assert i_i is not None
                          try:
                            sigmas.append(float(st[si_i]))
                          except:
                            proceed = False
                       flags_.append(st[t_i])
           if(not proceed):
              print >> err, "Could not get column data (hkl,fo,flags,sigma): ", hkl
           if(indices.size() == 0 or data.size() == 0 or len(flags_) == 0):
              proceed = False
           if(indices.size() > 0 and proceed):
              assert indices.size() == data.size()
              assert indices.size() == len(flags_)
              if(sigmas.size() > 0):
                 assert sigmas.size() == data.size()
              else:
                 sigmas = flex.double(data.size(), 1.0)
              flags, cif_selection, proceed = get_array_of_r_free_flags(
                                                         proceed = proceed,
                                                         flags   = flags_,
                                                         cs      = cs,
                                                         indices = indices,
                                                         pdb     = pdb,
                                                         err     = err)
     #
     # now we're ready to set up a miller array (and may be even output it)
     #
     if(proceed):
        if(f_i is not None): data_label = "FOBS"
        else:                data_label = "IOBS"
        proceed, file_name = \
                    compose_and_write_mtz(indices          = indices,
                                          data             = data,
                                          sigmas           = sigmas,
                                          flags            = flags,
                                          cif_selection    = cif_selection,
                                          pdb              = pdb,
                                          err              = err,
                                          proceed          = proceed,
                                          cs               = cs,
                                          observation_type = data_label)

        if(proceed):
           result, fmodel, format = one_run(pdb       = pdb,
                                            hkl       = file_name,
                                            error_out = err)
           os.system("rm -rf %s"%file_name)
           if(result is not None):
              proceed, file_name = \
                    compose_and_write_mtz(indices          = fmodel.f_obs.indices(),
                                          data             = fmodel.f_obs.data(),
                                          sigmas           = fmodel.f_obs.sigmas(),
                                          flags            = fmodel.r_free_flags.data(),
                                          cif_selection    = None,
                                          pdb              = pdb,
                                          err              = err,
                                          proceed          = proceed,
                                          cs               = cs,
                                          observation_type = "FOBS")
              ofo = open("_result_"+pdb[-8:-4],"w")
              print >> ofo, result
              ofo.flush()
              ofo.close()
              err.flush()
  err.close()
  os.system("rm -rf "+hkl[:-2])


if(__name__ == "__main__"):
   run(args = sys.argv[1:])
