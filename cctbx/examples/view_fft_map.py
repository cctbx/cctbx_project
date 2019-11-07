"""
Loads the file map_coeff.pickle (see random_f_calc.py) and displays
the FFT map based on these coefficients in PyMOL. Also computes
and displays a list of peaks in the map.

Usage:

Setup the cctbx environment (e.g. source setpaths.csh) and
launch PyMOL from the command line. Inside PyMOL enter:

run view_fft_map.py
show_fft()
"""
from __future__ import absolute_import, division, print_function
from six.moves import range
from six.moves import zip

print("Loading module:", __name__)

# cctbx imports
from cctbx import maptbx
from libtbx import easy_pickle

# PyMOL imports
from chempy.map import Map
from pymol import cmd
from pymol import cgo

def show_map(unit_cell, map_covering_unit_cell, label, level):
  map_grid = map_covering_unit_cell.focus()
  print("map_grid:", map_grid)
  ucell_params = unit_cell.parameters()
  first = [0,0,0]
  last = [map_grid[i] + 1 for i in range(3)]
  c_obj_map = maptbx.as_CObjectZYX(
    map_unit_cell=map_covering_unit_cell,
    first=first,
    last=last,
    apply_sigma_scaling=True)
  map=Map()
  map.from_c_object(c_obj_map,'CObjectZYXfloat',
                    ucell_params[0:3], ucell_params[3:6],
                    list(map_grid), first, last)
  cmd.load_map(map, label+"_cell")
  print("map loaded into PyMol")
  cmd.isomesh(label+"_con", label+"_cell", level) # create mesh
  cmd.color('gray', label+"_cell") # color wire frame
  cmd.set('auto_zoom', '0')    # disable zooming
  cmd.set('ortho', '1')        # orthoscopic projects
  cmd.enable(label+"_cell") # put box around map object
  cmd.color('cyan', label+"_con")   # color mesh

def show_peaks(unit_cell, clusters, radius=2.0):
  go = []
  go.extend([cgo.COLOR, 1, 0, 0,])
  height0 = None
  for site,height in zip(clusters.sites(), clusters.heights()):
    print("%8.5f %8.5f %8.5f" % site, height)
    if (height0 == None): height0 = height
    go.extend(  [cgo.SPHERE]
              + list(unit_cell.orthogonalize(site))
              + [radius*height/height0])
  cmd.load_cgo(go, "peaks")

def show_fft(file="map_coeff.pickle", map_level=3.0):
  cmd.delete("all")
  map_coeff = easy_pickle.load(file)
  map_coeff.show_summary()
  fft_map = map_coeff.fft_map(
    symmetry_flags=maptbx.use_space_group_symmetry)
  print("map gridding:", fft_map.n_real())
  show_map(fft_map.unit_cell(), fft_map.real_map(), "fft_map", map_level)
  clusters = fft_map.tags().peak_search(
    parameters=maptbx.peak_search_parameters(
      min_distance_sym_equiv=3.0,
      max_clusters=10),
    map=fft_map.real_map()).all()
  show_peaks(fft_map.unit_cell(), clusters)
  cmd.zoom('all', 15.0) # zoom with additional border of 15 Ang.
  print()

if (__name__ == "pymol"):
  cmd.extend("show_fft", show_fft)
