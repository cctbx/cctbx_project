#! /bin/csh -fe
set verbose
curl http://cci.lbl.gov/cctbx_build/results/current/cctbx_bundle.selfx -o cctbx_bundle.selfx
perl cctbx_bundle.selfx 0
rm cctbx_bundle.selfx
rm cctbx_install_script.csh
mv cctbx_sources sources
cd sources
git clone https://github.com/cctbx/cctbx_project.git
cctbx_project/libtbx/development/move_obsolete_from_bundle.csh
