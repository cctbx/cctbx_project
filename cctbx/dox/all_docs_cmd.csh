#! /bin/csh -f
set verbose
rm -rf ../../htdocs/current_cvs
mkdir -p ../../htdocs/current_cvs/python
doxygen
cd rst
docutils_cmd.csh
echo $cwd
cd ../../../htdocs/current_cvs/python
libtbx.help -w `libtbx.show_dist_paths libtbx`
libtbx.help -w `libtbx.show_dist_paths boost_adaptbx`
libtbx.help -w `libtbx.show_dist_paths scitbx`
libtbx.help -w `libtbx.show_dist_paths cctbx`
libtbx.help -w `libtbx.show_dist_paths iotbx`
libtbx.help -w `libtbx.show_dist_paths mmtbx`
