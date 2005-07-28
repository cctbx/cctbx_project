from iotbx import pdb
import iotbx.pdb.interpretation
import sys, math, time, os
from mmtbx.tls.tls import *
import mmtbx

def get_resolution(st):
  result = None
  q1 = (st.count("REMARK   2 ")==1 and st.count("RESOLUTION")>0 and
        st.count("ANGSTROM")>0)
  ch = st.split()
  if(q1):
     try:
       result = float(ch[3])
       assert result < 100.0 and result > 0.01
     except:
       pass
  return result

def extract_resolution(remark_2_records):
  res_counter = 0
  resolutions = []
  for record in remark_2_records:
      result = get_resolution(record)
      if(result is not None):
         resolutions.append(result)
         res_counter +=1
  if(res_counter == 1):
     return resolutions
  else:
     return None
