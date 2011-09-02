#! /bin/csh -fe
# cd <html-root>
# mkdir jv395; cd jv395
# unzip javaviewFull-v3.95.zip
set verbose
set cctbx="`libtbx.show_dist_paths cctbx`"
touch crmXXX
rm -rf crm???
libtbx.parallel_simple --dirs=crm --command="cctbx.python $cctbx/sgtbx/direct_space_asu/check_redundancies.py --strip_grid 24 "'$(MULTI:230-1)'
rm -rf asu_gallery
mkdir asu_gallery
mv crm???/asu_gallery/* asu_gallery
rm -rf crm???
libtbx.parallel_simple --dirs=crm --command="cctbx.python $cctbx/sgtbx/direct_space_asu/check_redundancies.py --strip_grid 24 --plane_group "'$(MULTI:17-1)'
mv crm???/asu_gallery/* asu_gallery
cctbx.python "$cctbx/web/asu_gallery/jv_asu.py" --server cci.lbl.gov
cctbx.python "$cctbx/web/asu_gallery/jv_asu.py" --server cci.lbl.gov --plane_group
