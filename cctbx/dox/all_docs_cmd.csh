#! /bin/csh -f
set verbose
cd `libtbx.show_repository_paths`
rm -rf ../htdocs/current
mkdir -p ../htdocs/current/python
libtbx.extract_code_from_txt libtbx/libtbx/phil/doc.txt
docutils.rst2html doc.rst > libtbx_phil.html
rm doc.rst
mv libtbx_phil.html libtbx_phil_examples.py ../htdocs/
doxygen cctbx/dox/Doxyfile
cd cctbx/dox/rst
docutils_cmd.csh
echo $cwd
cd ../../../../htdocs/current/python
libtbx.help -w `libtbx.show_dist_paths libtbx`
libtbx.help -w `libtbx.show_dist_paths boost_adaptbx`
libtbx.help -w `libtbx.show_dist_paths scitbx`
libtbx.help -w `libtbx.show_dist_paths cctbx`
libtbx.help -w `libtbx.show_dist_paths iotbx`
libtbx.help -w `libtbx.show_dist_paths mmtbx`
cd ../../siena2005
./update_generated_files.csh
cd ../sbgrid2008
./update_generated_files.csh
