#! /bin/csh -f
set echo
docutils.buildhtml -stg
mkdir -p ../../../htdocs/current_cvs
mv *.html ../../../htdocs/current_cvs
cp *.txt *.css *.png ../../../htdocs/current_cvs
cd ../../../htdocs/current_cvs
rm -f index.html
ln -s introduction.html index.html
