example_of_tls_parameters_in_remark_3 = """\
REMARK   3
REMARK   3  TLS DETAILS
REMARK   3   NUMBER OF TLS GROUPS  :    2
REMARK   3
REMARK   3   TLS GROUP :     1
REMARK   3    NUMBER OF COMPONENTS GROUP :    1
REMARK   3    COMPONENTS        C SSSEQI   TO  C SSSEQI
REMARK   3    RESIDUE RANGE :   A     1        A   348
REMARK   3    ORIGIN FOR THE GROUP (A):  40.4920  20.4900  77.3450
REMARK   3    T TENSOR
REMARK   3      T11:   0.0445 T22:   0.0094
REMARK   3      T33:   0.0759 T12:   0.0073
REMARK   3      T13:  -0.0019 T23:  -0.0026
REMARK   3    L TENSOR
REMARK   3      L11:   0.4084 L22:   0.4087
REMARK   3      L33:   0.4695 L12:   0.0456
REMARK   3      L13:   0.0685 L23:  -0.0388
REMARK   3    S TENSOR
REMARK   3      S11:  -0.0063 S12:   0.0233 S13:   0.0090
REMARK   3      S21:   0.0417 S22:  -0.0071 S23:   0.0132
REMARK   3      S31:  -0.0008 S32:  -0.0002 S33:   0.0134
REMARK   3
REMARK   3   TLS GROUP :     2
REMARK   3    NUMBER OF COMPONENTS GROUP :    1
REMARK   3    COMPONENTS        C SSSEQI   TO  C SSSEQI
REMARK   3    RESIDUE RANGE :   B XYZW1-       B XY348+
REMARK   3    ORIGIN FOR THE GROUP (A):  16.4920  -1.1810  59.3130
REMARK   3    T TENSOR
REMARK   3      T11:   0.0412 T22:   0.0118
REMARK   3      T33:   0.1002 T12:  -0.0082
REMARK   3      T13:  -0.0117 T23:  -0.0290
REMARK   3    L TENSOR
REMARK   3      L11:   0.3189 L22:   0.3849
REMARK   3      L33:   0.8975 L12:   0.0578
REMARK   3      L13:  -0.0398 L23:   0.1013
REMARK   3    S TENSOR
REMARK   3      S11:   0.0415 S12:   0.0183 S13:  -0.0365
REMARK   3      S21:   0.0087 S22:  -0.0145 S23:   0.0507
REMARK   3      S31:   0.1029 S32:  -0.0744 S33:  -0.0269
REMARK   3
"""

class tls(object):
   def __init__(self, T, L, S, origin, selection_string = None):
     self.T = T
     self.L = L
     self.S = S
     self.origin = origin
     self.selection_string = selection_string #XXX do this smarter

def extract_tls_parameters(remark_3_records, file_name = ""):
# T = (T11, T22, T33, T12, T13, T23)
# L = (L11, L22, L33, L12, L13, L23)
# S = (S11, S12, S13, S21, S22, S23, S31, S32, S33)
  T = []
  L = []
  S = []
  origin = []
  residue_range = []
  number_of_tls_groups = 0
  group_counter = 0
  record_start = False
  record_end = False
  one_tls_group_records = []
  all_tls_group_records = []
  tls_params = []
  for record in remark_3_records:
    assert record[0:10] == "REMARK   3"
    if(record.startswith("REMARK   3   TLS GROUP :")):
       group_number = None
       try: group_number = int(record.split()[5])
       except ValueError:
         print "Cannot extract TLS group number:\n  ", record, file_name
         return []
       record_start = True
    if(record.startswith("REMARK   3      S31:")):
       record_end = True
       record_start = False
       one_tls_group_records.append(record)
       all_tls_group_records.append(one_tls_group_records)
       one_tls_group_records = []
    else:
       record_end = False
    if(record_start and not record_end):
       one_tls_group_records.append(record)
  for one in all_tls_group_records:
    n_components = 0
    r_range = []
    origin = None
    T = None
    T11,T22,T33,T12,T13,T23 = None,None,None,None,None,None
    L11,L22,L33,L12,L13,L23 = None,None,None,None,None,None
    S11,S12,S13,S21,S22,S23,S31,S32,S33= None,None,None,None,None,None,None,None,None
    sel_str = None
    for rec in one:
      if(rec.startswith("REMARK   3    NUMBER OF COMPONENTS GROUP :")):
         n_components = None
         try: n_components = int(rec.split()[7])
         except ValueError:
           print "Cannot extract number of TLS components:\n  ",n_components,\
             file_name
           return []
      if(rec.startswith("REMARK   3    RESIDUE RANGE :")):
        ch1,res1,ch2,res2 = rec[32:33], rec[34:40], rec[47:48], rec[49:55]
        r_range.append([ch1,res1,ch2,res2])
      if(rec.startswith("REMARK   3    ORIGIN FOR THE GROUP (A):")):
         try: x = float(rec.split()[7])
         except ValueError:
           print "Cannot extract x of origin:\n  ", rec, file_name
           return []
         try: y = float(rec.split()[8])
         except ValueError:
           print "Cannot extract y of origin:\n  ", rec, file_name
           return []
         try: z = float(rec.split()[9])
         except ValueError:
           print "Cannot extract z of origin:\n  ", rec, file_name
           return []
         origin = [x,y,z]
      if(rec.startswith("REMARK   3      T11:")):
         assert [T11, T22, T33, T12, T13, T23].count(None) == 6
         try: T11 = float(rec.split()[3])
         except ValueError:
           print "Cannot extract T11:\n  ", rec, file_name
           return []
         try: T22 = float(rec.split()[5])
         except ValueError:
           print "Cannot extract T22:\n  ", rec, file_name
           return []
      if(rec.startswith("REMARK   3      T33:")):
         assert [T11, T22, T33, T12, T13, T23].count(None) == 4
         try: T33 = float(rec.split()[3])
         except ValueError:
           print "Cannot extract T33:\n  ", rec, file_name
           return []
         try: T12 = float(rec.split()[5])
         except ValueError:
           print "Cannot extract T12:\n  ", rec, file_name
           return []
      if(rec.startswith("REMARK   3      T13:")):
         assert [T11, T22, T33, T12, T13, T23].count(None) == 2
         try: T13 = float(rec.split()[3])
         except ValueError:
           print "Cannot extract T13:\n  ", rec, file_name
           return []
         try: T23 = float(rec.split()[5])
         except ValueError:
           print "Cannot extract T23:\n  ", rec, file_name
           return []
         T=[T11, T22, T33, T12, T13, T23]
      if(rec.startswith("REMARK   3      L11:")):
         assert [L11, L22, L33, L12, L13, L23].count(None) == 6
         try: L11 = float(rec.split()[3])
         except:
           try: L11 = float(rec[20:30])
           except ValueError:
             print "Cannot extract L11:\n  ", rec, file_name
             return []
         try: L22 = float(rec.split()[5])
         except:
           try: L22 = float(rec[34:44])
           except ValueError:
             print "Cannot extract L22:\n  ", rec, file_name
             return []
      if(rec.startswith("REMARK   3      L33:")):
         assert [L11, L22, L33, L12, L13, L23].count(None) == 4
         try: L33 = float(rec.split()[3])
         except:
           try: L33 = float(rec[20:30])
           except ValueError:
             print "Cannot extract L33:\n  ", rec, file_name
             return []
         try: L12 = float(rec.split()[5])
         except:
           try: L12 = float(rec[34:44])
           except ValueError:
             print "Cannot extract L12:\n  ", rec, file_name
             return []
      if(rec.startswith("REMARK   3      L13:")):
         assert [L11, L22, L33, L12, L13, L23].count(None) == 2
         try: L13 = float(rec.split()[3])
         except:
           try: L13 = float(rec[20:30])
           except ValueError:
             print "Cannot extract L13:\n  ", rec, file_name
             return []
         try: L23 = float(rec.split()[5])
         except:
           try: L23 = float(rec[34:44])
           except ValueError:
             print "Cannot extract L23:\n  ", rec, file_name
             return []
         L=[L11, L22, L33, L12, L13, L23]
      if(rec.startswith("REMARK   3      S11:")):
         assert [S11, S12, S13, S21, S22, S23, S31, S32, S33].count(None) == 9
         try: S11 = float(rec.split()[3])
         except ValueError:
           print "Cannot extract S11:\n  ", rec, file_name
           return []
         try: S12 = float(rec.split()[5])
         except ValueError:
           print "Cannot extract S12:\n  ", rec, file_name
           return []
         try: S13 = float(rec.split()[7])
         except ValueError:
           print "Cannot extract S13:\n  ", rec, file_name
           return []
      if(rec.startswith("REMARK   3      S21:")):
         assert [S11, S12, S13, S21, S22, S23, S31, S32, S33].count(None) == 6
         try: S21 = float(rec.split()[3])
         except ValueError:
           print "Cannot extract S21:\n  ", rec, file_name
           return []
         try: S22 = float(rec.split()[5])
         except ValueError:
           print "Cannot extract S22:\n  ", rec, file_name
           return []
         try: S23 = float(rec.split()[7])
         except ValueError:
           print "Cannot extract S23:\n  ", rec, file_name
           return []
      if(rec.startswith("REMARK   3      S31:")):
         assert [S11, S12, S13, S21, S22, S23, S31, S32, S33].count(None) == 3
         try: S31 = float(rec.split()[3])
         except ValueError:
           print "Cannot extract S31:\n  ", rec, file_name
           return []
         try: S32 = float(rec.split()[5])
         except ValueError:
           print "Cannot extract S32:\n  ", rec, file_name
           return []
         try: S33 = float(rec.split()[7])
         except:
           try:
             if(rec.split()[7].count("NULL")):
                S33 = - (S11 + S22)
           except ValueError:
             print "Cannot extract S33:\n  ", rec, file_name
             return []
         S=[S11, S12, S13, S21, S22, S23, S31, S32, S33]
    for rr in r_range:
      ch1 = rr[0].strip()
      ch2 = rr[2].strip()
      res1 = rr[1].strip()
      res2 = rr[3].strip()
      if(rr[0]==rr[2]):
        if(sel_str is None):
          if(ch1 == ""):
            sel_str = "(resid %s:%s)" % (res1,res2)
          else:
            sel_str = "(chain %s and resid %s:%s)" % (ch1,res1,res2)
        else:
          if(ch1 == ""):
            sel_str += " or (resid %s:%s)" % (res1,res2)
          else:
            sel_str += " or (chain %s and resid %s:%s)" % (ch1,res1,res2)
      elif(rr[0] != rr[2]):
        if(sel_str is None):
          if(ch1 == ""):
            sel_str = "(resid %s- ) or (resid -%s)" % (res1,res2)
          else:
            sel_str = "(chain %s and resid %s- ) or (chain %s and resid -%s)"%(
              ch1,res1,ch2,res2)
        else:
          if(ch1 == ""):
            sel_str += " or (resid %s- ) or (resid -%s)" % (res1,res2)
          else:
            sel_str += \
              " or (chain %s and resid %s- ) or (chain %s and resid -%s)" %(
              ch1,res1,ch2,res2)
    if(sel_str is not None or
       T.count(None) > 0 or
       L.count(None) > 0 or
       S.count(None) > 0):
      tls_params.append(tls(T=T,L=L,S=S,origin=origin,selection_string=sel_str))
  return tls_params


def get_program(st):
  st = st.upper()
  program = None
  programs = ["SHELX", "CNS", "PROLSQ", "FROG", "REFMAC", "XPLOR", "X-PLOR",
              "NUCLSQ", "MOPRO", "MOLLY", "PROFFT", "TNT", "CORELS", "GROMOS"
              "XTALVIEW", "RESTRAIN", "GPRLSA", "ARP", "CNX", "EREF",
              "NMREF", "GSAS", "BUSTER","SOLVE", "CCP4", "NUCLIN", "MAIN",
              "PHENIX", "PHENIX.REFINE"]
  diamond = "DIAMOND"
  if(st[0:6] == "REMARK"):
     st_split = st.split()
     ch_1 = st.count("REMARK   3") == 1
     ch_2 = st.count("PROGRAM") == 1
     ch_3 = st.count("REFINEMENT") == 1
     if(ch_1 and ch_2):
        for item in st_split:
            for potential_program in programs:
                if(potential_program.count(",")):
                   potential_program = potential_program[:-1]
                if(item.count(potential_program)):
                   program = potential_program
     if(program is None and ch_1):
        refinement_or_program = st.count("REFINEMENT")>=1 or \
                                st.count("PROGRAM")>=1 or \
                                st.count("PROCEDURE")>=1
        j_and_l = st.count("JACK")==1 or st.count("LEVITT")
        if(refinement_or_program and j_and_l): program = "JACK_AND_LEVITT"
        if((ch_2 or ch_3) and st.count(diamond) != 0): program = diamond
        k_and_h = st.count("KONNERT")==1 or st.count("HENDRICKSON")==1
        k_and_h_1 = st.count("KONNERT AND HENDRICKSON")==1
        if((refinement_or_program and k_and_h) or k_and_h_1):
           program = "KONNERT_AND_HENDRICKSON"
        sa = st.count("SIMULATED")==1 and st.count("ANNEALING")
        if(refinement_or_program and sa): program = "X-PLOR"
        if(refinement_or_program and st.count("CORELS")==1): program = "CORELS"
        jones_and_liljas = "JONES_AND_A.LILJAS"
        if(st.count(jones_and_liljas) or st.count("T.A.JONES, L.LILJAS")):
           program = jones_and_liljas
  if(program == "XPLOR"): program = "X-PLOR"
  return program

def format_name(program_names):
  try:
    new = program_names[0]
    first = program_names[0]
    for name in program_names:
      if(name != first):
         new +="/"+name
         first = name
  except: return program_names
  return new


def extract_program_name(remark_3_records):
  p_counter = 0
  program_names = []
  for record in remark_3_records:
      result = get_program(record)
      if(result is not None):
         program_names.append(result)
         p_counter +=1
  if(p_counter != 0):
     return format_name(program_names)
  else:
     return None
