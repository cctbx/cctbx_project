
from __future__ import absolute_import, division, print_function

from libtbx.utils import to_str

def load_refinement(
    pdb_file,
    map_file,
    work_dir=None,
    peaks_file=None,
    show_symmetry=True,
    have_anom_map=None,
    have_residual_map=None):
  import coot # import dependency
  import os.path
  import sys
  if (work_dir is not None):
    if os.path.isdir(work_dir):
      os.chdir(work_dir)
    else :
      for k, arg in enumerate(sys.argv[1:-1]):
        if (arg == "--script"):
          script_file = sys.argv[k+2]
          if os.path.isfile(script_file):
            script_dir = os.path.dirname(script_file)
            if (script_dir != ""):
              os.chdir(script_dir)
  read_pdb(to_str(pdb_file))
  if (peaks_file is not None) and (os.path.isfile(peaks_file)):
    imol = read_pdb(to_str(peaks_file))
    set_mol_displayed(imol, False)
  auto_read_make_and_draw_maps(to_str(map_file))
  if (have_anom_map):
    imol = make_and_draw_map(to_str(map_file),"ANOM","PHANOM","",0,0)
    set_contour_level_in_sigma(imol, 3.0)
    if (have_residual_map):
      set_map_colour(imol, 0.0, 1.0, 1.0)
      set_map_displayed(imol, False)
    else :
      set_map_colour(imol, 1.0, 1.0, 0.0)
  if (have_residual_map):
    imol = make_and_draw_map(to_str(map_file),"ANOMDIFF","PHANOMDIFF","",0,1)
    set_contour_level_in_sigma(imol, 3.0)
    set_map_colour(imol, 1.0, 1.0, 0.0)
