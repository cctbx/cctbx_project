#! /bin/csh -fe
set verbose
mkdir asu_gallery
cctbx.python "`libtbx.show_dist_paths cctbx`"/cctbx/sgtbx/direct_space_asu/check_redundancies.py --strip_grid 24 1-230
