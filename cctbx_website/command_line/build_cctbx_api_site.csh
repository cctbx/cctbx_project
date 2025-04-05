#!/bin/csh -f
echo "Building documentation API for cctbx_project"
setenv module_list  "fftw3tbx scitbx gltbx serialtbx chiltbx iota clipper_adaptbx iotbx simtbx cma_es kokkostbx smtbx cootbx libtbx sphinx crys3d mmtbx spotfinder boost cudatbx tbxx boost_adaptbx dox omptbx ucif cbflib_adaptbx dox.sphinx prime wxtbx cctbx fable qttbx xfel fast_linalg rstbx"
rm -fr working
mkdir working
cd working

foreach x ($module_list)
phenix.python ../run_pdoc_cctbx_api.py $x >& $x.log &
end
wait
echo ""
echo "Results by module:"
grep failed *.log | grep -v List

cd ..

echo ""
echo "Making directory cctbx_project_api with html"

if (-d cctbx_project_api) rm -rf cctbx_project_api
mkdir cctbx_project_api
cp cctbx_api_site_index.html cctbx_project_api/index.html
foreach x ($module_list)
  if (-d working/$x) mv working/$x cctbx_project_api/$x
end
tar czvf - cctbx_project_api > cctbx_project_api.tgz
ls -tlr cctbx_project_api
echo "Ready with cctbx_project_api in cctbx_project_api.tgz"
ls -tlr cctbx_project_api.tgz
echo "ALL DONE"
 
