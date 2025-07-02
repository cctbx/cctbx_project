from __future__ import absolute_import, division, print_function
from cctbx import maptbx
from scitbx.array_family import flex
from cctbx.eltbx import xray_scattering

import boost_adaptbx.boost.python as bp
ext = bp.import_ext("cctbx_maptbx_bcr_bcr_ext")

std_labels = xray_scattering.standard_labels_list()

def rfactor(a,b):
  n = flex.sum(flex.abs(a-b))
  d = flex.sum(flex.abs(a+b))
  return n/d*100*2

def run():
  o = maptbx.atom_curves(scattering_type="C", scattering_table="wk1995")
  b = o.bcr_approx(
    d_min       = 2.0,
    radius_max  = 12,
    radius_step = 0.01,
    kprot=2,
    mxp=1000,
    epsc=0.001,
    )
  r = rfactor(b.image_values, b.bcr_approx_values)
  print(r)
  for R, B, C in zip(b.R, b.B, b.C):
    print("%20.16f %20.16f %20.16f"%(R, B, C) )

if (__name__ == "__main__"):
  run()


#{
#  "2.0": {
#    "err": 0.0009001,
#    "R (M)": [...],
#    "B (N)": [...],
#    "C (K)": [...]
#  },
#  "2.5": {
#    "err": 0.0007166,
#    "R (M)": [...],
#    "B (N)": [...],
#    "C (K)": [...]
#  },
#  ...
#}
#
#
#import json
#import re
#
#text = """(your full text block here as a multi-line string)"""
#
#blocks = re.split(r'\bRes=', text)
#results = {}
#
#for block in blocks[1:]:  # Skip the first empty block before 'Res='
#    header, *rows = block.strip().split('\n')
#
#    # Extract Res and err
#    match = re.search(r'(\d+\.?\d*)\s+err=([\d.eE+-]+)', header)
#    if not match:
#        continue
#    res = float(match.group(1))
#    err = float(match.group(2))
#
#    # Skip the column header line
#    rows = rows[1:]
#
#    r_vals, b_vals, c_vals = [], [], []
#    for row in rows:
#        parts = row.strip().split()
#        if len(parts) != 3:
#            continue
#        r, b, c = map(float, parts)
#        r_vals.append(r)
#        b_vals.append(b)
#        c_vals.append(c)
#
#    results[str(res)] = {  # Use string key for JSON compliance
#        "err": err,
#        "R (M)": r_vals,
#        "B (N)": b_vals,
#        "C (K)": c_vals
#    }
#
## Write to JSON
#with open("data.json", "w") as f:
#    json.dump(results, f, indent=2)
