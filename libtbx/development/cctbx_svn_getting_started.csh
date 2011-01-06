#! /bin/csh -fe
set verbose
wget http://cci.lbl.gov/cctbx_build/results/current/cctbx_bundle.selfx
perl cctbx_bundle.selfx 0
rm cctbx_bundle.selfx
rm cctbx_install_script.csh
mv cctbx_sources sources
cd sources
svn co https://cctbx.svn.sourceforge.net/svnroot/cctbx/trunk cctbx_project
cctbx_project/libtbx/development/move_obsolete_from_bundle.csh
