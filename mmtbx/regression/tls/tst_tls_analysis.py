from __future__ import absolute_import, division, print_function
from mmtbx.tls import analysis
import os
import libtbx.load_env
from scitbx import matrix
from libtbx.test_utils import approx_equal

def extract(file_name):
  of = open(file_name, "r")
  start  = False
  start_T=False
  start_L=False
  start_S=False
  T=[]
  L=[]
  S=[]
  for l in of.readlines():
    l=l.strip()
    if(l.count("TLS BUILT FROM THE BASE ELEMENTS")):
      start=True
    if(start):
      if(l.startswith("matrix T")):
        start_T=True
      if(l.startswith("matrix L")):
        start_L=True
      if(l.startswith("matrix S")):
        start_S=True
    if(start_T and len(l.split())==3):
      if(len(T)<3): T.append(l)
      else:         start_T=False
    if(start_L and len(l.split())==3):
      if(len(L)<3): L.append(l)
      else:         start_L=False
    if(start_S and len(l.split())==3):
      if(len(S)<3): S.append(l)
      else:         start_S=False
  of.close()
  def convert(m):
    return matrix.sqr([
      float(m[0].split()[0]), float(m[0].split()[1]), float(m[0].split()[2]),
      float(m[1].split()[0]), float(m[1].split()[1]), float(m[1].split()[2]),
      float(m[2].split()[0]), float(m[2].split()[1]), float(m[2].split()[2])
    ])
  T = convert(T)
  L = convert(L)
  S = convert(S)
  return T, L, S

files = [
  "dec04_test000.mes",
  "dec04_test001.mes",
  "dec04_test002.mes",
  "dec04_test003.mes",
  "dec04_test004.mes",
  "dec04_test005.mes",
  "dec04_test012.mes",
  "dec04_test013.mes",
  "dec04_test014.mes",
  "dec04_test015.mes",
  "dec04_test021.mes",
  "dec04_test022.mes",
  "dec04_test023.mes",
  "dec04_test024.mes",
  "dec04_test033.mes",
  "dec04_test034.mes",
  "dec04_test035.mes",
  "dec04_test044.mes",
  "dec04_test113.mes",
  "dec04_test114.mes",
  "dec04_test115.mes",
  "dec04_test144.mes",
  "dec04_test145.mes",
  "dec04_test146.mes",
  "dec04_test147.mes",
  "dec04_test148.mes",
  "dec04_test149.mes",
  "dec04_test159.mes",

  "dec04_test168.mes",
  "dec04_test169.mes",
  "dec04_test215.mes"
]

def run():
  for fn in files:
    try:
      fn_ = libtbx.env.find_in_repositories(
        relative_path="mmtbx/regression/tls/data_tls_analysis/%s"%fn,
        test=os.path.isfile)
      of = open("phenix_"+fn, "w")
      T, L, S = extract(fn_)
      # test python implementation
      r = analysis.run(T=T, L=L, S=S, log=of, implementation='python')
      # test c++ implementation
      r2 = analysis.run(T=T, L=L, S=S, log=of, implementation='c++')
      # test equivalence
      assert approx_equal(r.result.dx, r2.result.dx, 1e-4)
      assert approx_equal(r.result.dy, r2.result.dy, 1e-4)
      assert approx_equal(r.result.dz, r2.result.dz, 1e-4)
      assert approx_equal(r.result.l_x, r2.result.l_x, 1e-4)
      assert approx_equal(r.result.l_y, r2.result.l_y, 1e-4)
      assert approx_equal(r.result.l_z, r2.result.l_z, 1e-4)
      assert approx_equal(r.result.sx, r2.result.sx, 1e-4)
      assert approx_equal(r.result.sy, r2.result.sy, 1e-4)
      assert approx_equal(r.result.sz, r2.result.sz, 1e-4)
      assert approx_equal(r.result.tx, r2.result.tx, 1e-4)
      assert approx_equal(r.result.ty, r2.result.ty, 1e-4)
      assert approx_equal(r.result.tz, r2.result.tz, 1e-4)
      assert approx_equal(r.result.v_x, r2.result.v_x, 1e-4)
      assert approx_equal(r.result.v_x_M, r2.result.v_x_M, 1e-4)
      assert approx_equal(r.result.v_y, r2.result.v_y, 1e-4)
      assert approx_equal(r.result.v_y_M, r2.result.v_y_M, 1e-4)
      assert approx_equal(r.result.v_z, r2.result.v_z, 1e-4)
      assert approx_equal(r.result.v_z_M, r2.result.v_z_M, 1e-4)
      assert approx_equal(r.result.w_L_lx, r2.result.w_L_lx, 1e-4)
      assert approx_equal(r.result.w_L_ly, r2.result.w_L_ly, 1e-4)
      assert approx_equal(r.result.w_L_lz, r2.result.w_L_lz, 1e-4)
      assert approx_equal(r.result.w_M_lx, r2.result.w_M_lx, 1e-4)
      assert approx_equal(r.result.w_M_ly, r2.result.w_M_ly, 1e-4)
      assert approx_equal(r.result.w_M_lz, r2.result.w_M_lz, 1e-4)

      of.close()
      print(fn, "OK")
    except Exception as e:
      print(fn, str(e))


if (__name__ == "__main__"):
  run()
