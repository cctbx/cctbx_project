from iotbx import pdb
import iotbx.pdb.interpretation
import sys, math, time, os
from mmtbx.tls.tls import *
import mmtbx

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
  assert len(T) == len(L) == len(S) == len(origin) == len(residue_range)
  return mmtbx.tls.tls.tls_parameters(
                        T = T,
                        L = L,
                        S = S,
                        origin = origin,
                        number_of_tls_groups = number_of_tls_groups,
                        residue_range = residue_range)
