#! /bin/csh -fe
set verbose
set html_dir=/net/boa/srv/html/cci/gltbx
cp README.txt $html_dir
cd $html_dir
docutils.rst2html -stg README.txt README.html
