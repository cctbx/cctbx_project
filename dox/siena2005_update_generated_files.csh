#! /bin/csh -fe
source "`libtbx.show_build_path`/setpaths_all.csh"
set PHENIX_REGRESSION_DIR="`libtbx.find_in_repositories phenix_regression`"
set verbose
#
cp -r $CCTBX_DIST/examples/unit_cell_refinement.txt .
docutils.rst2html -stg unit_cell_refinement.txt > unit_cell_refinement.html
cctbx.python $CCTBX_DIST/examples/unit_cell_refinement.py > unit_cell_refinement.out
#
cp -r $IOTBX_DIST/examples/direct_methods_light.txt .
docutils.rst2html -stg direct_methods_light.txt > direct_methods_light.html
cp -r "$PHENIX_REGRESSION_DIR/misc/"{vj1132Isup2.hkl,vj1132sup1.cif} .
iotbx.python $IOTBX_DIST/examples/direct_methods_light.py vj1132Isup2.hkl vj1132sup1.cif > direct_methods_light.out
