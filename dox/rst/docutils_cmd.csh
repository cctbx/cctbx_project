#! /bin/csh -f
set echo
docutils.buildhtml -stg
mkdir -p ../../../../htdocs/current
mv *.html ../../../../htdocs/current
cp *.txt *.css *.png ../../../../htdocs/current
cd ../../../../htdocs/current
rm -f index.html
ln -s introduction.html index.html
