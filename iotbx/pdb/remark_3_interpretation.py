import os
from libtbx import group_args

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
   def __init__(self, t, l, s, origin, selection_string = None):
     self.t = t
     self.l = l
     self.s = s
     self.origin = origin
     self.selection_string = selection_string #XXX do this smarter

def has_tls(remark_3_records):
  result = False
  n_t = 0
  n_l = 0
  n_s = 0
  for line in remark_3_records:
    if(line.count("T11") or line.count("T 11")): n_t += 1
    if(line.count("T22") or line.count("T 22")): n_t += 1
    if(line.count("T33") or line.count("T 33")): n_t += 1
    if(line.count("T12") or line.count("T 12")): n_t += 1
    if(line.count("T13") or line.count("T 13")): n_t += 1
    if(line.count("T23") or line.count("T 23")): n_t += 1
    #
    if(line.count("L11") or line.count("L 11")): n_l += 1
    if(line.count("L22") or line.count("L 22")): n_l += 1
    if(line.count("L33") or line.count("L 33")): n_l += 1
    if(line.count("L12") or line.count("L 12")): n_l += 1
    if(line.count("L13") or line.count("L 13")): n_l += 1
    if(line.count("L23") or line.count("L 23")): n_l += 1
    #
    if(line.count("S11") or line.count("S 11")): n_s += 1
    if(line.count("S21") or line.count("S 21")): n_s += 1
    if(line.count("S31") or line.count("S 31")): n_s += 1
    if(line.count("S12") or line.count("S 12")): n_s += 1
    if(line.count("S22") or line.count("S 22")): n_s += 1
    if(line.count("S32") or line.count("S 32")): n_s += 1
    if(line.count("S13") or line.count("S 13")): n_s += 1
    if(line.count("S23") or line.count("S 23")): n_s += 1
    if(line.count("S33") or line.count("S 33")): n_s += 1
  if(n_t > 3 and n_l > 3 and n_s > 3): result = True
  return result

def prepocess_line(line):
  l0 = line.split()
  new_elements = []
  for l_ in l0:
    if(l_.isalpha() or l_.isdigit()): new_elements.append(l_)
    else:
      try:
        val = float(l_)
        new_elements.append(l_)
      except:
        tmp = ""
        for i, c in enumerate(l_):
          if(i == 0): tmp+=c
          elif(c in ["+","-"]):
            if(not l_[i-1].isalpha()):
              new_elements.append(tmp)
              tmp = c
            else: tmp+=c
          else: tmp+=c
        new_elements.append(tmp)
  return " ".join(new_elements)

class extract_tls_parameters(object):

   def __init__(self, remark_3_records, chain_ids, pdb_hierarchy, file_name = ""):
     # T = (T11, T22, T33, T12, T13, T23)
     # L = (L11, L22, L33, L12, L13, L23)
     # S = (S11, S12, S13, S21, S22, S23, S31, S32, S33)
     self.pdb_hierarchy = pdb_hierarchy
     self.remark_3_records = remark_3_records
     self.chain_ids = chain_ids
     self.file_name = os.path.basename(file_name)
     self.tls_present = has_tls(remark_3_records = remark_3_records)
     self.tls_params = []
     self.error_string = None
     if(self.tls_present):
       self.extract()
     if(self.error_string is not None):
       self.tls_params = []

   def format_err(self, msg, rec=""):
     x = self.file_name.strip()+" : "
     y = " : "+rec.strip()
     if(len(y.strip())==1): y = ""
     if(len(x.strip())==1): x = ""
     self.error_string = x+msg.strip()+y

   def extract(self):
     T = []
     L = []
     S = []
     origin = []
     residue_range = []
     group_counter = 0
     record_start = False
     record_end = False
     one_tls_group_records = []
     all_tls_group_records = []
     for record in self.remark_3_records:
       assert record[0:10] == "REMARK   3"
       if(record.startswith("REMARK   3   TLS GROUP :")):
          group_number = None
          try: group_number = int(record.split()[5])
          except ValueError:
            self.format_err(msg="Cannot extract TLS group number:", rec=record)
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
       T11,T22,T33,T12,T13,T23 = [None]*6
       L11,L22,L33,L12,L13,L23 = [None]*6
       S11,S12,S13,S21,S22,S23,S31,S32,S33 = [None]*9
       sel_str = None
       for i_seq, rec in enumerate(one):
         if(rec.startswith("REMARK   3    NUMBER OF COMPONENTS GROUP :")):
            n_components = None
            try: n_components = int(rec.split()[7])
            except ValueError:
              self.format_err(msg="Cannot extract number of TLS components.", rec=rec)
              return []
         if(rec.startswith("REMARK   3    RESIDUE RANGE :")):
           if(len(rec.split())==9):
             rec_s = rec.split()
             ch1,res1,ch2,res2 = rec_s[5],rec_s[6],rec_s[7],rec_s[8]
           else:
             ch1,res1,ch2,res2 = rec[31:33], rec[34:40], rec[47:48], rec[49:55]
           r_range.append([ch1,res1,ch2,res2])
         # PHENIX selection
         if(rec.startswith("REMARK   3    SELECTION:")):
           sel_str = rec[rec.index(":")+1:]
           if(sel_str.strip().upper() in ["NONE","NULL"]):
             self.format_err(msg="Bad TLS selection string1.", rec = rec)
             return []
           i = i_seq+1
           while ( one[i].startswith("REMARK   3             :") or
                   one[i].startswith("REMARK   3              ") ):
             if(one[i].count("ORIGIN")): break
             sel_str += " "+one[i][24:]
             i += 1
             if(sel_str.strip().upper() in ["NONE","NULL"]):
               self.format_err(msg="Bad TLS selection string2.", rec = rec)
               return []
           sel_str = " ".join(sel_str.split())
           ##
           sel_str_spl = sel_str.split()
           new_str = ""
           for ie, e in enumerate(sel_str_spl):
             if(ie==0): new_str = e
             else:
               if(new_str[len(new_str)-1]==":"):
                 if(e[0].isdigit()):
                   new_str += e
                 else:
                   new_str = new_str + " " + e
               elif(new_str[len(new_str)-1].isdigit()):
                 if(e[0]==":"):
                   new_str += e
                 else:
                   new_str = new_str + " " + e
               else:
                 new_str = new_str + " " + e
           sel_str = new_str
           ##
         #
         if(rec.startswith("REMARK   3    ORIGIN FOR THE GROUP (A):")):
            rec = prepocess_line(rec)
            try: x = float(rec.split()[7])
            except ValueError:
              self.format_err(msg="Cannot extract origin.", rec=rec)
              return []
            try: y = float(rec.split()[8])
            except ValueError:
              self.format_err(msg="Cannot extract origin.", rec=rec)
              return []
            try: z = float(rec.split()[9])
            except ValueError:
              self.format_err(msg="Cannot extract origin.", rec=rec)
              return []
            origin = [x,y,z]
         if(rec.startswith("REMARK   3      T11:")):
            assert [T11, T22, T33, T12, T13, T23].count(None) == 6
            try: T11 = float(rec.split()[3])
            except ValueError:
              self.format_err(msg="Cannot extract T.", rec=rec)
              return []
            try: T22 = float(rec.split()[5])
            except ValueError:
              self.format_err(msg="Cannot extract T.", rec=rec)
              return []
         if(rec.startswith("REMARK   3      T33:")):
            assert [T11, T22, T33, T12, T13, T23].count(None) == 4
            try: T33 = float(rec.split()[3])
            except ValueError:
              self.format_err(msg="Cannot extract T.", rec=rec)
              return []
            try: T12 = float(rec.split()[5])
            except ValueError:
              self.format_err(msg="Cannot extract T.", rec=rec)
              return []
         if(rec.startswith("REMARK   3      T13:")):
            assert [T11, T22, T33, T12, T13, T23].count(None) == 2
            try: T13 = float(rec.split()[3])
            except ValueError:
              self.format_err(msg="Cannot extract T.", rec=rec)
              return []
            try: T23 = float(rec.split()[5])
            except ValueError:
              self.format_err(msg="Cannot extract T.", rec=rec)
              return []
            T=[T11, T22, T33, T12, T13, T23]
         if(rec.startswith("REMARK   3      L11:")):
            assert [L11, L22, L33, L12, L13, L23].count(None) == 6
            try: L11 = float(rec.split()[3])
            except:
              try: L11 = float(rec[20:30])
              except ValueError:
                self.format_err(msg="Cannot extract L.", rec=rec)
                return []
            try: L22 = float(rec.split()[5])
            except:
              try: L22 = float(rec[34:44])
              except ValueError:
                self.format_err(msg="Cannot extract L.", rec=rec)
                return []
         if(rec.startswith("REMARK   3      L33:")):
            assert [L11, L22, L33, L12, L13, L23].count(None) == 4
            try: L33 = float(rec.split()[3])
            except:
              try: L33 = float(rec[20:30])
              except ValueError:
                self.format_err(msg="Cannot extract L.", rec=rec)
                return []
            try: L12 = float(rec.split()[5])
            except:
              try: L12 = float(rec[34:44])
              except ValueError:
                self.format_err(msg="Cannot extract L.", rec=rec)
                return []
         if(rec.startswith("REMARK   3      L13:")):
            assert [L11, L22, L33, L12, L13, L23].count(None) == 2
            try: L13 = float(rec.split()[3])
            except:
              try: L13 = float(rec[20:30])
              except ValueError:
                self.format_err(msg="Cannot extract L.", rec=rec)
                return []
            try: L23 = float(rec.split()[5])
            except:
              try: L23 = float(rec[34:44])
              except ValueError:
                self.format_err(msg="Cannot extract L.", rec=rec)
                return []
            L=[L11, L22, L33, L12, L13, L23]
         if(rec.startswith("REMARK   3      S11:")):
            assert [S11, S12, S13, S21, S22, S23, S31, S32, S33].count(None) == 9
            try: S11 = float(rec.split()[3])
            except ValueError:
              self.format_err(msg="Cannot extract S.", rec=rec)
              return []
            try: S12 = float(rec.split()[5])
            except ValueError:
              self.format_err(msg="Cannot extract S.", rec=rec)
              return []
            try: S13 = float(rec.split()[7])
            except ValueError:
              self.format_err(msg="Cannot extract S.", rec=rec)
              return []
         if(rec.startswith("REMARK   3      S21:")):
            assert [S11, S12, S13, S21, S22, S23, S31, S32, S33].count(None) == 6
            try: S21 = float(rec.split()[3])
            except ValueError:
              self.format_err(msg="Cannot extract S.", rec=rec)
              return []
            try: S22 = float(rec.split()[5])
            except ValueError:
              self.format_err(msg="Cannot extract S.", rec=rec)
              return []
            try: S23 = float(rec.split()[7])
            except ValueError:
              self.format_err(msg="Cannot extract S.", rec=rec)
              return []
         if(rec.startswith("REMARK   3      S31:")):
            assert [S11, S12, S13, S21, S22, S23, S31, S32, S33].count(None) == 3
            try: S31 = float(rec.split()[3])
            except ValueError:
              self.format_err(msg="Cannot extract S.", rec=rec)
              return []
            try: S32 = float(rec.split()[5])
            except ValueError:
              self.format_err(msg="Cannot extract S.", rec=rec)
              return []
            try: S33 = float(rec.split()[7])
            except:
              try:
                if(rec.split()[7].count("NULL")):
                   S33 = - (S11 + S22)
              except ValueError:
                self.format_err(msg="Cannot extract S.", rec=rec)
                return []
            S=[S11, S12, S13, S21, S22, S23, S31, S32, S33]
       if(len(r_range) > 0):
         assert sel_str is None
         for rr in r_range:
           ch1 = rr[0].strip()
           ch2 = rr[2].strip()
           res1 = rr[1].strip()
           res2 = rr[3].strip()
           if(len(res1)>0 and len(res2)>0):
             try: res1_ = int(res1)
             except ValueError:
               self.format_err(msg="Bad TLS selection string3.", rec = "".join(rr))
               return []
             try: res2_ = int(res2)
             except ValueError:
               self.format_err(msg="Bad TLS selection string4.", rec = "".join(rr))
               return []
             if(res1_ > res2_ and ch1 == ch2):
               self.format_err(msg="Bad TLS selection: start index > end index.")
               return []
           tmp = refmac_range_to_phenix_string_selection(
             pdb_hierarchy = self.pdb_hierarchy, chain_start = ch1,
             resseq_start = res1, chain_end = ch2, resseq_end = res2)
           if(sel_str is None): sel_str = tmp
           else: sel_str = sel_str +" or %s"%tmp
       if(sel_str is not None):
         if(sel_str.count("(") != sel_str.count(")")):
           self.format_err(msg="Bad TLS selection: missing ) or (", rec = sel_str)
           return []
         if(sel_str.upper() in ["NONE","NULL"]):
           sel_str = None
           self.format_err(msg="Bad TLS selection string5.", rec = sel_str)
           return []
       if(sel_str is not None or
          T.count(None) > 0 or
          L.count(None) > 0 or
          S.count(None) > 0):
         # Compatibility with previous versions of PHENIX: replace "-" with ":"
         new_c = ""
         for i,c in enumerate(sel_str):
           try: cl = sel_str[i-1]
           except: cl = c
           if(c=="-" and cl.isdigit()): c = ":"
           new_c += c
         sel_str = new_c
         #
         self.tls_params.append(tls(
           t=T,l=L,s=S,origin=origin,selection_string=sel_str))
     if(self.tls_present and len(self.tls_params)==0):
       self.format_err(msg="TLS present but cannot be extracted.")
       return []
     #
     if(len(self.tls_params)>0):
       chain_ids_new_u = [c.upper() for c in self.chain_ids]
       chain_ids_new_l = [c.lower() for c in self.chain_ids]
       if(chain_ids_new_u!=self.chain_ids and chain_ids_new_l!=self.chain_ids):
         self.format_err(msg="Mixed chain ids detected.")
         return []
     #
     if(len(all_tls_group_records) != len(self.tls_params)):
       self.format_err(msg="TLS present but cannot be extracted.")
       return []
     #
     for t in self.tls_params:
       if(not(len(t.t) == len(t.l))):
         self.format_err(msg="TLS present but cannot be extracted.")
         return []
       if(len(t.t) == len(t.l) and len(t.l)>0 and
          (len(t.s)!=9 or len(t.origin)!=3)):
         self.format_err(msg="TLS present but cannot be extracted.")
         return []
     #
     return self.tls_params


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

def refmac_range_to_phenix_string_selection(pdb_hierarchy, chain_start,
      resseq_start, chain_end, resseq_end):
  chain_start  =  chain_start.strip()
  resseq_start =  resseq_start.strip()
  chain_end    =  chain_end.strip()
  resseq_end   =  resseq_end.strip()
  sel_str = None
  if(chain_start == chain_end):
    if(chain_start != ""):
      if([resseq_start,resseq_end].count("")==0):
        sel_str = "(chain %s and resseq %s:%s)"%(chain_start,resseq_start,resseq_end)
      else:
        sel_str = "(chain %s)"%(chain_start)
    else:
      if([resseq_start,resseq_end].count("")==0):
        sel_str = "(resseq %s:%s)"%(resseq_start,resseq_end)
  else:
    sel_str1 = None
    if(chain_start != ""):
      if(resseq_start != ""):
        sel_str1 = "(chain %s and resseq %s:)"%(chain_start,resseq_start)
      else:
        sel_str1 = "(chain %s)"%(chain_start)
    else:
      if(resseq_start != ""):
        sel_str1 = "(resseq %s:)"%(resseq_start)
    sel_str2 = None
    if(chain_end != ""):
      if(resseq_end != ""):
        sel_str2 = "(chain %s and resseq :%s)"%(chain_end,resseq_end)
      else:
        sel_str2 = "(chain %s)"%(chain_end)
    else:
      if(resseq_end != ""):
        sel_str2 = "(resseq :%s)"%(resseq_end)
    if([sel_str1,sel_str2].count(None)==0):
      sel_str = "%s or %s"%(sel_str1, sel_str2)
    elif(sel_str1 is not None): sel_str = sel_str1
    elif(sel_str2 is not None): sel_str = sel_str2
    if(sel_str is not None):
      models = pdb_hierarchy.models()
      if(len(models)>1): return None # XXX one model only
      start_collecting = False
      chain_ids = []
      for chain in models[0].chains():
        chain_id = chain.id.strip()
        if(chain_id == chain_end): break
        if(chain_id == chain_start):
          start_collecting = True
          continue
        if(start_collecting): chain_ids.append(chain_id)
      if(len(chain_ids)>0):
        for chain_id in chain_ids:
          sel_str = sel_str+" or chain %s"%chain_id
  if(sel_str is not None):
    return "%s"%sel_str
  else: return None


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

def extract_f_model_core_constants(remark_3_records):
  # XXX phenix.refine specific
  k_sol            = None
  b_sol            = None
  b_cart           = [None]*6
  twin_fraction    = None
  twin_law         = None
  r_solv           = None
  r_shrink         = None
  grid_step_factor = None
  b11,b22,b33,b12,b13,b23 = [None,None,None,None,None,None]
  r_work           = None
  r_free           = None
  for l in remark_3_records:
    l = l.strip()
    ls= l.split()
    try:
      if(l.count("REMARK   3   K_SOL              :")): k_sol = float(ls[4])
      if(l.count("REMARK   3   B_SOL              :")): b_sol = float(ls[4])
      if(l.count("REMARK   3    B11")): b11 = float(ls[len(ls)-1])
      if(l.count("REMARK   3    B22")): b22 = float(ls[len(ls)-1])
      if(l.count("REMARK   3    B33")): b33 = float(ls[len(ls)-1])
      if(l.count("REMARK   3    B12")): b12 = float(ls[len(ls)-1])
      if(l.count("REMARK   3    B13")): b13 = float(ls[len(ls)-1])
      if(l.count("REMARK   3    B23")): b23 = float(ls[len(ls)-1])
      if(l.count("REMARK   3   FRACTION:")): twin_fraction = float(ls[3])
      if(l.count("REMARK   3   OPERATOR:")): twin_law      = ls[3]
      if(l.count("REMARK   3   SOLVENT RADIUS     :")): r_solv = float(ls[5])
      if(l.count("REMARK   3   SHRINKAGE RADIUS   :")): r_shrink = float(ls[5])
      if(l.count("REMARK   3   GRID STEP FACTOR   :")): grid_step_factor = float(ls[6])
      if(l.count("REMARK   3   R VALUE            (WORKING SET) :")): r_work = float(ls[7])
      if(l.count("REMARK   3   FREE R VALUE                     :")): r_free = float(ls[6])
    except: pass
  if([b11,b22,b33,b12,b13,b23].count(None)==0): b_cart=[b11,b22,b33,b12,b13,b23]
  return group_args(
    k_sol            = k_sol,
    b_sol            = b_sol,
    b_cart           = b_cart,
    twin_fraction    = twin_fraction,
    twin_law         = twin_law,
    r_solv           = r_solv,
    r_shrink         = r_shrink,
    grid_step_factor = grid_step_factor,
    r_work           = r_work,
    r_free           = r_free)
