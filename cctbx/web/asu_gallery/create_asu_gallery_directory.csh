#! /bin/csh -fe
# cd <html-root>
# mkdir jv395; cd jv395
# unzip javaviewFull-v3.95.zip
set verbose
mkdir asu_gallery
cctbx.python "`libtbx.show_dist_paths cctbx`"/sgtbx/direct_space_asu/check_redundancies.py --strip_grid 24 1-230
cctbx.python "`libtbx.show_dist_paths cctbx`"/web/asu_gallery/jv_asu.py --server cci.lbl.gov
