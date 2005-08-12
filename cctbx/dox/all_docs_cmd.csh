#! /bin/csh -f
rm -rf ../../htdocs/current_cvs
mkdir -p ../../htdocs/current_cvs
doxygen
cd rst
docutils_cmd.csh
