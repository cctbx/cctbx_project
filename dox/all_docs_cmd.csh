#! /bin/csh -f
set verbose
cd "`libtbx.show_dist_paths libtbx`"
cd ..
cd ..
set root="`pwd`"
if (! -d "$root/cctbx_project") then
  echo "ERROR: unexpected directory structure: no $root/cctbx_project directory."
  exit 1
endif
if (! -d "$root/htdocs") then
  echo "ERROR: unexpected directory structure: no $root/htdocs directory."
  exit 1
endif
#
if (1) then
  cd "$root/htdocs"
  libtbx.extract_code_from_txt "$root/cctbx_project/libtbx/phil/doc.txt"
  docutils.rst2html doc.rst > libtbx_phil.html
  rm doc.rst
endif
#
if (1) then
  cd "$root/htdocs"
  rm -rf current/c_plus_plus
  mkdir -p current
  cd "$root/cctbx_project"
  doxygen dox/Doxyfile
endif
#
if (1) then
  cd "$root/cctbx_project/dox/rst"
  docutils.buildhtml -stg
  mv *.html "$root/htdocs/current"
  cp *.txt *.css *.png "$root/htdocs/current"
  cd "$root/htdocs/current"
  rm -f index.html
  ln -s introduction.html index.html
  cd "$root/htdocs"
  rm -f index.html
  ln -s current/versions.html index.html
endif
#
if (1) then
  mkdir -p "$root/htdocs/current"
  cd "$root/htdocs/current"
  rm -rf python
  mkdir python
  cd python
  libtbx.help -w "$root/cctbx_project/boost_adaptbx"
  libtbx.help -w "$root/cctbx_project"
endif
#
if (1) then
  cd "$root/htdocs"
  rm -rf siena2005
  mkdir siena2005
  cd siena2005
  cp -p "$root/cctbx_project/dox/siena2005/"* .
  cp -p "$root/cctbx_project/dox/rst/default.css" .
  "$root/cctbx_project/dox/siena2005_update_generated_files.csh"
endif
#
if (1) then
  cd "$root/htdocs"
  rm -rf sbgrid2008
  mkdir sbgrid2008
  cd sbgrid2008
  cp -r "$root/cctbx_project/dox/rst/default.css" .
  "$root/cctbx_project/dox/sbgrid2008_update_generated_files.csh"
endif
