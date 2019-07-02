from __future__ import absolute_import, division, print_function
import os,sys
from cctbx.examples.merging import test_levenberg_sparse as test
import libtbx.load_env

# test script assumes that you get the data files directly from the author (NKS) and
# install them in the directory "xscale_reserve" at the same dir-level as cctbx_project

if __name__=="__main__":
  modules_dist = os.path.abspath(os.path.join(libtbx.env.dist_path("cctbx"),"../.."))
  datadir = os.path.join(modules_dist,"xscale_reserve") # Get files directly from author, NKS
  plot_flag=False
  esd_plot_flag=False
  for N in [25, 200, 300, 400, 500, 800, 1000, 2000, 5000]:
   for trans in [1.0, 0.1, 0.01, 0.001, 0.0001, 0.00001]:
    test.execute_case(datadir, n_frame=N, transmittance=trans, apply_noise=True,
      plot=plot_flag, esd_plot = esd_plot_flag)
    print("OK")
    #raw_input("OK")
    sys.stdout.flush()
