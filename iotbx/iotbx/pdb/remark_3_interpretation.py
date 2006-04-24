from iotbx import pdb
import iotbx.pdb.interpretation
import sys, math, time, os

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
REMARK   3    RESIDUE RANGE :   B     1        B   348
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

def extract_tls_parameters(remark_3_records):
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
       try: group_number = int(record.split()[5])
       except ValueError: print "Cannot extract TLS group number."
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
    sel_string = None
    for rec in one:
      if(rec.startswith("REMARK   3    NUMBER OF COMPONENTS GROUP :")):
         try: n_components = int(rec.split()[7])
         except ValueError: print "Cannot extract number of TLS components."
      if(rec.startswith("REMARK   3    RESIDUE RANGE :")):
         if(len(rec.split()) == 9):
            rr = rec.split()
            chain_first = rr[5]
            try: residue_number_first = int(rr[6])
            except ValueError: print "Cannot extract first residue number in residue range."
            chain_second = rr[7]
            #if(chain_first != chain_second):
            #   raise RuntimeError("chain_first != chain_second: %s %s" % (chain_first,chain_second))
            try: residue_number_second = int(rr[8])
            except ValueError: print "Cannot extract second residue number in residue range."
            r_range.append([chain_first,residue_number_first,
                            chain_second,residue_number_second])
         elif(len(rec.split()) == 7):
            rr = rec.split()
            chain_first = " "
            try: residue_number_first = int(rr[5])
            except ValueError: print "Cannot extract first residue number in residue range."
            chain_second = " "
            try: residue_number_second = int(rr[6])
            except ValueError: print "Cannot extract second residue number in residue range."
            r_range.append([chain_first,residue_number_first,
                            chain_second,residue_number_second])
      if(rec.startswith("REMARK   3    ORIGIN FOR THE GROUP (A):")):
         try: x = float(rec.split()[7])
         except ValueError: print "Cannot extract x of origin."
         try: y = float(rec.split()[8])
         except ValueError: print "Cannot extract y of origin."
         try: z = float(rec.split()[9])
         except ValueError: print "Cannot extract z of origin."
         origin = [x,y,z]
      if(rec.startswith("REMARK   3      T11:")):
         assert [T11, T22, T33, T12, T13, T23].count(None) == 6
         try: T11 = float(rec.split()[3])
         except ValueError: print "Cannot extract T11."
         try: T22 = float(rec.split()[5])
         except ValueError: print "Cannot extract T22."
      if(rec.startswith("REMARK   3      T33:")):
         assert [T11, T22, T33, T12, T13, T23].count(None) == 4
         try: T33 = float(rec.split()[3])
         except ValueError: print "Cannot extract T33."
         try: T12 = float(rec.split()[5])
         except ValueError: print "Cannot extract T12."
      if(rec.startswith("REMARK   3      T13:")):
         assert [T11, T22, T33, T12, T13, T23].count(None) == 2
         try: T13 = float(rec.split()[3])
         except ValueError: print "Cannot extract T13."
         try: T23 = float(rec.split()[5])
         except ValueError: print "Cannot extract T23."
         T=[T11, T22, T33, T12, T13, T23]
      if(rec.startswith("REMARK   3      L11:")):
         assert [L11, L22, L33, L12, L13, L23].count(None) == 6
         try: L11 = float(rec.split()[3])
         except:
           try: L11 = float(rec[20:30])
           except ValueError: print "Cannot extract L11."
         try: L22 = float(rec.split()[5])
         except:
           try: L22 = float(rec[34:44])
           except ValueError: print "Cannot extract L22."
      if(rec.startswith("REMARK   3      L33:")):
         assert [L11, L22, L33, L12, L13, L23].count(None) == 4
         try: L33 = float(rec.split()[3])
         except:
           try: L33 = float(rec[20:30])
           except ValueError: print "Cannot extract L33."
         try: L12 = float(rec.split()[5])
         except:
           try: L12 = float(rec[34:44])
           except ValueError: print "Cannot extract L12."
      if(rec.startswith("REMARK   3      L13:")):
         assert [L11, L22, L33, L12, L13, L23].count(None) == 2
         try: L13 = float(rec.split()[3])
         except:
           try: L13 = float(rec[20:30])
           except ValueError: print "Cannot extract L13."
         try: L23 = float(rec.split()[5])
         except:
           try: L23 = float(rec[34:44])
           except ValueError: print "Cannot extract L23."
         L=[L11, L22, L33, L12, L13, L23]
      if(rec.startswith("REMARK   3      S11:")):
         assert [S11, S12, S13, S21, S22, S23, S31, S32, S33].count(None) == 9
         try: S11 = float(rec.split()[3])
         except ValueError: print "Cannot extract S11."
         try: S12 = float(rec.split()[5])
         except ValueError: print "Cannot extract S12."
         try: S13 = float(rec.split()[7])
         except ValueError: print "Cannot extract S13."
      if(rec.startswith("REMARK   3      S21:")):
         assert [S11, S12, S13, S21, S22, S23, S31, S32, S33].count(None) == 6
         try: S21 = float(rec.split()[3])
         except ValueError: print "Cannot extract S21."
         try: S22 = float(rec.split()[5])
         except ValueError: print "Cannot extract S22."
         try: S23 = float(rec.split()[7])
         except ValueError: print "Cannot extract S23."
      if(rec.startswith("REMARK   3      S31:")):
         assert [S11, S12, S13, S21, S22, S23, S31, S32, S33].count(None) == 3
         try: S31 = float(rec.split()[3])
         except ValueError: print "Cannot extract S31."
         try: S32 = float(rec.split()[5])
         except ValueError: print "Cannot extract S32."
         try: S33 = float(rec.split()[7])
         except:
           try:
             if(rec.split()[7].count("NULL")):
                S33 = - (S11 + S22)
           except ValueError: print "Cannot extract S33."
         S=[S11, S12, S13, S21, S22, S23, S31, S32, S33]
    for rr in r_range:
      if(rr[0] != " " and rr[0]==rr[2]):
         if(sel_string is None):
           sel_string = "(chain "+rr[0] + " and resid "+str(rr[1]) + ":"+str(rr[3])+")"
         else:
           sel_string += " or (chain "+rr[0] + " and resid "+str(rr[1]) + ":"+str(rr[3])+")"
      elif(rr[0] != rr[2]):
         if(sel_string is None):
           sel_string = "(chain "+rr[0] + " and resid "+str(rr[1]) + "-"+")" + \
                        " or (chain "+rr[2] + " and resid -"+str(rr[3]) +")"
         else:
           sel_string = " or (chain "+rr[0] + " and resid "+str(rr[1]) + "-"+")" + \
                        " or (chain "+rr[2] + " and resid -"+str(rr[3]) +")"
      else:
         if(sel_string is None):
           sel_string = "(resid "+str(rr[1]) + ":"+str(rr[3])+")"
         else:
           sel_string += " or (resid "+str(rr[1]) + ":"+str(rr[3])+")"

    tls_params.append(tls(T=T,L=L,S=S,origin=origin,selection_string=sel_string))
  return tls_params


def XXXextract_tls_parameters(remark_3_records):
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
  for record in remark_3_records:
    assert record[0:10] == "REMARK   3"
    srecord = record.split()
    lsrecord = len(srecord)
#=====> extract number of TLS groups
    if(record.count("NUMBER OF TLS GROUPS")*record.count(":") == 1):
      try: number_of_tls_groups = int(srecord[7])
      except:
        try: number_of_tls_groups = int(srecord[6])
        except:
          try: number_of_tls_groups = int(srecord[lsrecord-1])
          except:
            if(srecord[lsrecord-1][0] == ":"):
              try: number_of_tls_groups = int(srecord[lsrecord-1][1:])
              except ValueError: return None
            elif(srecord[lsrecord-1].count("S:") == 1):
              start_position = srecord[lsrecord-1].index(":")+1
              try: number_of_tls_groups = int(srecord[lsrecord-1][start_position:])
              except ValueError: return None
            else:
              return None
#=====> extract currents TLS group number
    if(record.count("TLS GROUP")*record.count(":") == 1 and
       record.count("NUMBER OF TLS GROUPS") == 0):
      try:
        group_number = int(srecord[5])
        group_counter += 1
      except:
        try:
          group_number = int(srecord[4])
          group_counter += 1
        except:
          try:
            group_number = int(srecord[lsrecord-1])
            group_counter += 1
          except:
            if(srecord[lsrecord-1][0] == ":"):
              try:
                group_number = int(srecord[lsrecord-1][1:])
                group_counter += 1
              except ValueError: return None
            elif(srecord[lsrecord-1].count("P:") == 1):
              start_position = srecord[lsrecord-1].index(":")+1
              try:
                group_number = int(srecord[lsrecord-1][start_position:])
                group_counter += 1
              except ValueError: return None
            else:
              return None
#=====> extract current residue range
    if(record.count("RESIDUE RANGE ")*record.count(":") == 1):
      rec = list(srecord[5:9])
      try:
        rec[1]
      except IndexError: return None
      try:
        int(rec[1])
      except ValueError: return None
      try:
        rec[3]
      except IndexError: return None
      try:
        int(rec[3])
      except ValueError: return None
      residue_range.append(list(srecord[5:9]))
    elif(record.count("RESIDUE RANGE:") == 1):
      rec = list(srecord[4:8])
      try:
        rec[1]
      except IndexError: return None
      try:
        int(rec[1])
      except ValueError: return None
      try:
        rec[3]
      except IndexError: return None
      try:
        int(rec[3])
      except ValueError: return None
      residue_range.append(list(srecord[4:8]))
#=====> extract current origin
    if(record.count("ORIGIN FOR THE GROUP") == 1):
      try: x = float(srecord[lsrecord-3])
      except ValueError: return None
      try: y = float(srecord[lsrecord-2])
      except ValueError: return None
      try: z = float(srecord[lsrecord-1])
      except ValueError: return None
      origin.append([x,y,z])
#=====> extract current T
    check = record.count("T11:")*record.count("T22:")
    if(check == 1 and lsrecord == 6):
      T11=None
      T22=None
      T33=None
      T12=None
      T13=None
      T23=None
      try: T11 = float(srecord[lsrecord-3])
      except ValueError: return None
      try: T22 = float(srecord[lsrecord-1])
      except ValueError: return None
    check = record.count("T33:")*record.count("T12:")
    if(check == 1 and lsrecord == 6):
      assert T11 is not None and T22 is not None
      assert T33 is None and T12 is None and T13 is None and T23 is None
      try: T33 = float(srecord[lsrecord-3])
      except ValueError: return None
      try: T12 = float(srecord[lsrecord-1])
      except ValueError: return None
    check = record.count("T13:")*record.count("T23:")
    if(check == 1 and lsrecord == 6):
      assert T11 is not None and T22 is not None
      assert T33 is not None and T12 is not None
      assert T13 is None and T23 is None
      try: T13 = float(srecord[lsrecord-3])
      except ValueError: return None
      try: T23 = float(srecord[lsrecord-1])
      except ValueError: return None
      T.append([T11, T22, T33, T12, T13, T23])
#=====> extract current L
    check = record.count("L11:")*record.count("L22:")
    if(check == 1 and lsrecord == 6):
      L11=None
      L22=None
      L33=None
      L12=None
      L13=None
      L23=None
      try: L11 = float(srecord[lsrecord-3])
      except ValueError: return None
      try: L22 = float(srecord[lsrecord-1])
      except ValueError: return None
    check = record.count("L33:")*record.count("L12:")
    if(check == 1 and lsrecord == 6):
      assert L11 is not None and L22 is not None
      assert L33 is None and L12 is None and L13 is None and L23 is None
      try: L33 = float(srecord[lsrecord-3])
      except ValueError: return None
      try: L12 = float(srecord[lsrecord-1])
      except ValueError: return None
    check = record.count("L13:")*record.count("L23:")
    if(check == 1 and lsrecord == 6):
      assert L11 is not None and L22 is not None
      assert L33 is not None and L12 is not None
      assert L13 is None and L23 is None
      try: L13 = float(srecord[lsrecord-3])
      except ValueError: return None
      try: L23 = float(srecord[lsrecord-1])
      except ValueError: return None
      L.append([L11, L22, L33, L12, L13, L23])
#=====> extract current S
    check = record.count("S11:")*record.count("S12:")*record.count("S13:")
    if(check == 1 and lsrecord == 8):
      S11=None
      S12=None
      S13=None
      S21=None
      S22=None
      S23=None
      S31=None
      S32=None
      S33=None
      try: S11 = float(srecord[lsrecord-5])
      except ValueError: return None
      try: S12 = float(srecord[lsrecord-3])
      except ValueError: return None
      try: S13 = float(srecord[lsrecord-1])
      except ValueError: return None
    check = record.count("S21:")*record.count("S22:")*record.count("S23:")
    if(check == 1 and lsrecord == 8):
      assert S11 is not None and S12 is not None and S13 is not None
      assert S21 is None and S22 is None and S23 is None
      assert S31 is None and S32 is None and S33 is None
      try: S21 = float(srecord[lsrecord-5])
      except ValueError: return None
      try: S22 = float(srecord[lsrecord-3])
      except ValueError: return None
      try: S23 = float(srecord[lsrecord-1])
      except ValueError: return None
    check = record.count("S31:")*record.count("S32:")*record.count("S33:")
    if(check == 1 and lsrecord == 8):
      assert S11 is not None and S12 is not None and S13 is not None
      assert S21 is not None and S22 is not None and S23 is not None
      assert S31 is None and S32 is None and S33 is None
      try: S31 = float(srecord[lsrecord-5])
      except ValueError: return None
      try: S32 = float(srecord[lsrecord-3])
      except ValueError: return None
      try: S33 = float(srecord[lsrecord-1])
      except ValueError: return None
      S.append([S11, S12, S13, S21, S22, S23, S31, S32, S33])

  assert group_counter == number_of_tls_groups == len(T)
  #assert len(T) == len(L) == len(S) == len(origin) == len(residue_range)
  #XXX self.pdb_residue_range: do this smarter
  tls_params = []
  for t,l,s,o,rr in zip(T,L,S,origin,residue_range):
      tls_params.append(tls(T=t,L=l,S=s,origin=o,pdb_residue_range=rr))
  return tls_params

def get_program(st):
  program = None
  programs = ["SHELX", "CNS", "PROLSQ", "FROG", "REFMAC", "XPLOR", "X-PLOR",
              "NUCLSQ", "MOPRO", "MOLLY", "PROFFT", "TNT", "CORELS", "GROMOS"
              "XTALVIEW", "RESTRAIN", "GPRLSA", "ARP", "CNX", "EREF",
              "NMREF", "GSAS", "BUSTER","SOLVE", "CCP4", "NUCLIN", "MAIN"]
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
        if(refinement_or_program and j_and_l): program = "JACK AND LEVITT"
        if((ch_2 or ch_3) and st.count(diamond) != 0): program = diamond
        k_and_h = st.count("KONNERT")==1 or st.count("HENDRICKSON")==1
        k_and_h_1 = st.count("KONNERT AND HENDRICKSON")==1
        if((refinement_or_program and k_and_h) or k_and_h_1):
           program = "KONNERT AND HENDRICKSON"
        sa = st.count("SIMULATED")==1 and st.count("ANNEALING")
        if(refinement_or_program and sa): program = "X-PLOR"
        if(refinement_or_program and st.count("CORELS")==1): program = "CORELS"
        jones_and_liljas = "JONES AND A. LILJAS"
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
