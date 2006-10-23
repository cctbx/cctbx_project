import sys, os, time
from cctbx.array_family import flex
from iotbx import crystal_symmetry_from_any
from cctbx import miller

models = "/net/cci-filer1/vol1/tmp/ftp.rcsb.org/uncompressed/pdb/"
data_files = "/net/cci-filer1/vol1/tmp/ftp.rcsb.org/all/structure_factors/"

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

tw_flags = ["1","1.","1.0","1.00","1.000","1.0000","-1","-1.","-1.0","-1.00",
            "-1.000","0","0.0","0.","0.00","0.000","0.0000","x","o","f","-",
            ">","<","h","l"]
cif_flags = ["x","o","f","-",">","<","h","l"]

ccp4_range = range(-99,-1)+ range(2,100)

def run():
  counter = 0
  nselected = 0
  hkl_file_names = os.listdir(data_files)
  pdb_file_names = os.listdir(models)
  ofo = open("_file_name_list","w")
  print "\n Number of data files  = ", len(hkl_file_names)
  print "\n Number of model files = ", len(pdb_file_names)
  #for hkl in ["r1a7qsf.ent.Z"]:
  for hkl in hkl_file_names:
      proceed = True
      data    = flex.double()
      sigmas  = flex.double()
      indices = flex.miller_index()
      flags_  = []
      flags = flex.bool()
      #
      # getting keys
      #
      if(hkl.count(".Z") == 1):
         counter += 1
         os.system("cp "+ data_files + hkl + " .")
         os.system("gunzip " + hkl)
         ifo = open(hkl[:-2], "r")
         lines = ifo.readlines()
         ifo.close()
         os.system("rm -rf *.ent*")
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
               keys = []
               proceed = False
               print "Multiple keys: ", keys, counter, hkl
               break
         if(proceed and len(keys) <= 4):
            print "No keys: ", keys, counter, hkl
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
                   print "\n Unknown key= ", key,  keys_updated, hkl, counter
            if(n_t > 1 and proceed):
               print "\n Multiple CV sets: ", keys_updated, hkl, counter, n_t
               proceed = False
            if(n_t == 0 and proceed):
               proceed = False
               print "\n No CV sets: ", keys_updated, hkl, counter, n_t
            if([f_i,i_i].count(None) == 2 or n_f > 1 or n_i > 1 and proceed):
               print "\n Multiple/no data sets: ", keys_updated,hkl,counter,n_f,n_i
               proceed = False
            if([h_i,k_i,l_i].count(None) > 0 or n_h > 1 or n_k > 1 or n_l > 1 and proceed):
               print "\n Multiple/no index sets: ",\
                      keys_updated,hkl,counter,n_h,n_k,n_l,"=",h_i,k_i,l_i
               proceed = False
         #
         # Hopefully, we know the numbers of colums with data. Try get it now.
         # Here we know the data is unique.
            if(proceed):
               for st in lines:
                  if(st == ''): break
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
               if(indices.size() == 0 or data.size() == 0 or len(flags_) == 0):
                  proceed = False
               if(indices.size() > 0 and proceed):
                  assert indices.size() == data.size()
                  assert indices.size() == len(flags_)
                  if(sigmas.size() > 0):
                     assert sigmas.size() == data.size()
                  else:
                     sigmas = flex.double(data.size(), 1.0)
                  flag_values = []
                  guess_ccp4_flags = 0
                  has_zero = False
                  guess_cif = 0
                  cif_selection = None
                  for item in flags_:
                      if(item not in flag_values):
                         flag_values.append(item)
                      try:
                         v = int(float(item))
                         if(v not in [-1,0,1]):
                            if(v in ccp4_range): guess_ccp4_flags += 1
                         if(v == 0): has_zero = True
                      except: continue
                      if(item in cif_flags): guess_cif += 1
                  if(guess_ccp4_flags > 0 and guess_cif > 0): proceed = False
                  if(proceed):
                     if(len(flag_values) == 2 and not guess_ccp4_flags and not guess_cif):
                        for item in flags_:
                            if(item == flag_values[0]): flags.append(True)
                            if(item == flag_values[1]): flags.append(False)
                     elif(len(flag_values) > 2 and guess_ccp4_flags and has_zero):
                        # we believe that "0" is "True" (test reflection)
                        for item in flags_:
                            if(int(float(item)) == 0): flags.append(True)
                            else:                      flags.append(False)
                     elif(guess_cif and "f" in flags_ and "o" in flags_):
                        cif_selection = flex.bool()
                        for item in flags_:
                            if(item == "f"):
                               flags.append(True)
                               cif_selection.append(True)
                            elif(item == "o"):
                               flags.append(False)
                               cif_selection.append(True)
                            else:
                               cif_selection.append(False)
                        if(cif_selection.size() != indices.size()):
                           proceed = False
                           print "\n Wrong flag type and number: ",hkl,counter,flag_values
                     else:
                        proceed = False
                        print "\n Wrong flag type: ",hkl,counter,flag_values
                     if(flags.count(True) > flags.count(False)):
                        flags = ~flags
                     if(flags.size() > 0):
                        n_trues = flags.count(True)
                        p_trues = n_trues*100./flags.size()
                        if(p_trues > 15.0 or p_trues < 2.0 or n_trues < 100):
                           print "\n Wrong flag amount = ",\
                                           hkl,counter,n_trues,p_trues,flag_values
                           proceed = False
         #
         # now we're ready to set up a miller array (and may be even output it)
         #
         if(proceed):
            nselected += 1
            write_mtz = True
            assert indices.size() == flags.size()
            assert sigmas.size() == data.size()
            assert data.size() == flags.size()
            assert flags.count(True) < flags.count(False)
            if(cif_selection is not None):
               indices = indices.select(cif_selection)
               sigmas  = sigmas.select(cif_selection)
               data    = data.select(cif_selection)
            assert indices.size() == flags.size()
            assert sigmas.size() == data.size()
            assert data.size() == flags.size()
            assert flags.count(True) < flags.count(False)
            for pdb in pdb_file_names:
                if(pdb.count(hkl[1:-8]) > 0):
                   break
            if(pdb.count(hkl[1:-8]) == 0):
               print "\n No PDB file found: ", pdb, hkl, nselected, counter
               write_mtz = False
            try:
              cs = crystal_symmetry_from_any.extract_from(models+pdb)
            except:
              write_mtz = False
              print "\n Cannot get cryst.symm.: ", pdb, hkl, nselected, counter
            try:
              miller_array = miller.set(
                               crystal_symmetry = cs,
                               indices          = indices,
                               anomalous_flag   = False).array(data   = data,
                                                               sigmas = sigmas)
              if(f_i is not None):
                 miller_array.set_observation_type_xray_amplitude()
                 data_label = "FOBS"
              else:
                 miller_array.set_observation_type_xray_intensity()
                 data_label = "IOBS"
            except:
              write_mtz = False
              print "\n Cannot setup miller arr.: ", pdb,hkl,nselected,counter
            if(write_mtz):
               try:
                 mtz_dataset = miller_array.as_mtz_dataset(
                                                column_root_label = data_label)
               except:
                 proceed = False
                 print "\n Cannot write MTZ file= ", pdb,hkl,nselected,counter
               if(proceed):
                  mtz_dataset.add_miller_array(
                          miller_array      = miller_array.array(data = flags),
                          column_root_label = "FLAGS")
                  mtz_object = mtz_dataset.mtz_object()
                  mtz_object.write(file_name = hkl[1:5]+".mtz")
                  d_max_min = miller_array.d_max_min()
                  size = miller_array.data().size()
                  print >> ofo, "%10s %10s %8.3f %8.3f %s %8.4f %10d" % (pdb,hkl[1:5]+".mtz",
                     d_max_min[0],d_max_min[1],miller_array.observation_type(),
                     flags.count(True)*100./flags.size(),size)
         os.system("rm -rf *ent*")

if (__name__ == "__main__"):
  t1 = time.time()
  run()
  t2 = time.time()
  print "\n Total time = ", t2 - t1
