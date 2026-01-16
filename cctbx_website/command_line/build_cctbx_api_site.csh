#!/bin/csh -f

# Supply either no arguments (runs cctbx API)
#   or NPROC, PREFIX, INDEX_FILE, extra module list

setenv module_list  "fftw3tbx scitbx gltbx serialtbx chiltbx iota clipper_adaptbx iotbx simtbx cma_es kokkostbx smtbx cootbx libtbx crys3d mmtbx spotfinder boost cudatbx tbxx boost_adaptbx dox omptbx ucif cbflib_adaptbx prime wxtbx cctbx fable qttbx xfel fast_linalg rstbx"

setenv base `libtbx.find_in_repositories cctbx_website`/command_line

if ($#argv > 0)then
  setenv NPROC $argv[1]
else
  setenv NPROC 4
endif

if ($#argv > 1) then
  setenv PREFIX $argv[2]
else
  setenv PREFIX cctbx_project
endif

if ($#argv > 2) then
  setenv INDEX_FILE $argv[3]
else
  setenv INDEX_FILE $base/cctbx_api_site_index.html
endif

if ($#argv > 3) then
    setenv module_list "$module_list $argv[4-$#argv]"
endif

echo "Building documentation API for $PREFIX using $NPROC processors"
echo "Base index.html will be from $INDEX_FILE"

echo "Module list: $module_list"

phenix.python -m pdoc --help > /dev/null
if ($status) then
  echo ""
  echo "This script needs pdoc3  nltk beautifulsoup4 . Install with:"
  echo ""
  echo "phenix.python -m pip install pdoc3 nltk beautifulsoup4"
  echo""
  exit 1
endif

# Work in a working directory that will be deleted at the end
if (-d working) then
  mv working working.old
endif

mkdir working
cd working

# Edit a few special files that cannot be interpreted as is by pdoc3

setenv cctbx_project `libtbx.find_in_repositories cctbx_project`
setenv files_to_edit "`libtbx.find_in_repositories iotbx`/pdb/hierarchy.py `libtbx.find_in_repositories scitbx`/array_family/flex.py"

@ i = 0
foreach x ($module_list)
@ i = $i + 1
#  Get HTML from python file, converting comments at top of functions and
#   classes to docstrings if no doc strings are present
phenix.python $base/run_pdoc_cctbx_api.py $x \
  convert_comments_to_docstring \
  files_to_edit_for_boost="${files_to_edit}" >& $x.log &
if ($i >= $NPROC) then
 @ i = 0
 wait
endif

end

wait

# Add the base html index.html
echo "Editing html files to simplify and add a base link"
cp $INDEX_FILE index.html

# Edit all the files to simplify and add a base link
foreach f (index.html */*.html */*/*.html */*/*/*.html */*/*/*/*.html)
  phenix.python $base/edit_html.py $f index_files $PREFIX >& /dev/null
end
echo ""
echo "Results by module:"
grep failed *.log | grep -v List

cd ..

setenv API_DIR ${PREFIX}_api
echo ""
echo "Making directory $API_DIR with html"

if (-d $API_DIR) rm -rf $API_DIR
mkdir $API_DIR
mv working/index.html $API_DIR/index.html
foreach x ($module_list)
  if (-d working/$x) mv working/$x $API_DIR/$x
end

# Add an index in $API_DIR/index_files
phenix.python $cctbx_project/libtbx/word_index_generator.py $API_DIR $API_DIR/index_files "$PREFIX API" > api_index.log

echo "Removing temporary files..."
if (-d working) rm -rf working
if (-f api_index.log ) rm -f api_index.log

ls -tlr $API_DIR
echo "Ready with $API_DIR"
echo "ALL DONE"

finish:
exit 0
