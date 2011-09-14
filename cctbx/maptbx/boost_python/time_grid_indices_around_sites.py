from __future__ import division
import iotbx.pdb
from cctbx import maptbx
from cctbx.array_family import flex
from scitbx import fftpack
from scitbx import matrix
from libtbx.utils import time_log
import sys

def run(args):
  assert len(args) == 1
  timer = time_log("pdb.input").start()
  pdb_inp = iotbx.pdb.input(file_name=args[0])
  print "number of pdb atoms:", pdb_inp.atoms().size()
  print timer.log()
  crystal_symmetry = pdb_inp.crystal_symmetry()
  assert crystal_symmetry is not None
  crystal_symmetry.show_summary()
  assert crystal_symmetry.unit_cell() is not None
  assert crystal_symmetry.space_group_info() is not None
  sites_cart = pdb_inp.atoms().extract_xyz()
  site_radii = flex.double(sites_cart.size(), 2.5)
  crystal_gridding = maptbx.crystal_gridding(
    unit_cell=crystal_symmetry.unit_cell(),
    d_min=2,
    resolution_factor=1/3)
  fft = fftpack.real_to_complex_3d(crystal_gridding.n_real())
  print "n_real:", fft.n_real()
  print "m_real:", fft.m_real()
  timer = time_log("grid_indices_around_sites").start()
  grid_indices = maptbx.grid_indices_around_sites(
    unit_cell=crystal_symmetry.unit_cell(),
    fft_n_real=fft.n_real(),
    fft_m_real=fft.m_real(),
    sites_cart=sites_cart,
    site_radii=site_radii)
  print "grid_indices.size():", grid_indices.size()
  print timer.log()
  print "grid fraction:", \
    grid_indices.size() / matrix.col(fft.n_real()).product()

if (__name__ == "__main__"):
  run(sys.argv[1:])
