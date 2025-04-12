#!/bin/csh -f
echo "Building documentation API for cctbx_project"
echo "WARNING: This version temporarily edits files in cctbx_project directory"
echo "Do not do anything in cctbx_project while this is running"

if (! -d $PHENIX/modules)then
  echo "This script needs PHENIX to be defined"
  goto finish
endif


setenv base $PHENIX/modules/cctbx_project/cctbx_website/command_line
setenv cctbx_project $PHENIX/modules/cctbx_project

setenv files_to_edit "iotbx/pdb/hierarchy.py scitbx/array_family/flex.py"

# Save original version of files to edit
foreach f ($files_to_edit)
  if (! -f $cctbx_project/$f)then
    echo "The file $cctbx_project/$f is missing"
    goto finish
  endif
  cp -p $cctbx_project/$f $cctbx_project/$f.original_version
  # Edit this file
  phenix.python $base/edit_for_boost.py $cctbx_project/$f
end

setenv module_list  "fftw3tbx scitbx gltbx serialtbx chiltbx iota clipper_adaptbx iotbx simtbx cma_es kokkostbx smtbx cootbx libtbx crys3d mmtbx spotfinder boost cudatbx tbxx boost_adaptbx dox omptbx ucif cbflib_adaptbx prime wxtbx cctbx fable qttbx xfel fast_linalg rstbx"
rm -fr working
mkdir working
cd working

foreach x ($module_list)
phenix.python $base/run_pdoc_cctbx_api.py $x >& $x.log &
end

echo ""
echo "WARNING: This version temporarily edits files in cctbx_project directory"
echo "Do not do anything in cctbx_project while this is running"
echo ""
wait

#Restore original files
foreach f ($files_to_edit)
  mv $cctbx_project/$f.original_version $cctbx_project/$f
enda

echo ""
echo "Original files in cctbx_project restored"
echo ""

# Add the base html index.html
echo "Editing html files to simplify and add a base link"

cp $base/cctbx_api_site_index.html index.html
# Edit all the files to simplify and add a base link
foreach f (index.html */*.html */*/*.html */*/*/*.html */*/*/*/*.html)
  phenix.python $base/edit_html.py $f index_files &
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
mv working/index.html cctbx_project_api/index.html
foreach x ($module_list)
  if (-d working/$x) mv working/$x cctbx_project_api/$x
end

# Add an index in cctbx_project_api/index_files
phenix.python $PHENIX/modules/cctbx_project/libtbx/word_index_generator.py cctbx_project_api cctbx_project_api/index_files "CCTBX" > api_index.log

echo "Packaging up files..."
tar czf - cctbx_project_api > cctbx_project_api.tgz
ls -tlr cctbx_project_api
echo "Ready with cctbx_project_api in cctbx_project_api.tgz"
ls -tlr cctbx_project_api.tgz
echo "ALL DONE"

finish:
