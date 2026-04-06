from __future__ import absolute_import, division, print_function
import libtbx
import gzip
import os, math
import json
from scitbx.array_family import flex

import boost_adaptbx.boost.python as bp
ext = bp.import_ext("cctbx_maptbx_bcr_bcr_ext")

def load_table(element=None, table=None, file_name=None):
  if file_name is None:
    element = element.strip().upper()
    path=libtbx.env.find_in_repositories("cctbx/maptbx/bcr/tables")
    file_name = "%s/%s_%s.json.gz"%(path, element, table)
  assert os.path.isfile(file_name)
  with gzip.open(file_name, "rt", encoding="utf-8") as f:
    return json.load(f)

def scatterers(xray_structure,
               resolution  = None,
               resolutions = None,
               RadFact     = 2.0,
               RadAdd      = 0.5):
  table = xray_structure.get_scattering_table()
  assert table in ["electron", "wk1995"]
  assert [resolution, resolutions].count(None)==1
  unit_cell = xray_structure.unit_cell()
  if resolutions is None:
    resolutions = [resolution,] * xray_structure.scatterers().size()
  RadMu   = RadFact + RadAdd
  arrays = {}
  element_types = list(
    xray_structure.scattering_type_registry().type_count_dict().keys())
  for e in element_types:
    d = load_table(element=e, table=table)
    arrays[e] = d
  ScaleB = 1.0 / (8.0 * math.pi**2)
  kscale = math.pi**1.5
  result = []
  shown = []
  for r, scatterer in zip(resolutions, xray_structure.scatterers()):
    e = scatterer.scattering_type.strip().upper()
    entry = arrays[e]
    keys = [float(x) for x in entry.keys()]
    key = str(min(keys, key=lambda x: abs(x - r)))
    vals = entry[key]
    R = flex.double(vals['R'])
    B = flex.double(vals['B'])
    C = flex.double(vals['C'])
    sel = R < (r*RadMu)
    R = R.select(sel)
    B = B.select(sel)
    C = C.select(sel)
    #if show_BCR and not e in shown:
    #  shown.append(e)
    #  print("    %s: R B C"%e)
    #  for r,b,c in zip(R,B,C):
    #    print("%15.9f %15.9f %15.9f"%(r,b,c))

    mu    = R
    nu    = B * ScaleB
    kappa = C
    musq  = mu * mu
    kappi = kappa/kscale

    bcr_scatterer = ext.bcr_scatterer(
      scatterer = scatterer,
      radius    = r*RadFact, # atomic radius = atomic_resolution * RadFact
      resolution=r,
      mu        = mu,
      kappa     = kappa,
      nu        = nu,
      musq      = musq,
      kappi     = kappi)
    result.append(bcr_scatterer)
  return result
